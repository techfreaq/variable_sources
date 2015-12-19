import os
import cPickle as pickle
from sklearn import cross_validation, ensemble
import sklearn
import pandas
import scipy
import numpy as N
import matplotlib.pyplot as plt

def dump_csc_ptf_features():
    """Write CSC and PTF variability features to a pickle file"""
    db = DB(os.environ['LSD_DB'])

    SIGMA_TO_95pctCL = 1.95996

    def max_radii(*args):
        return N.max(N.array(zip(args)),axis=0).flatten()

    with open('csc_feature_names.txt','r') as f:
        csc_keys = f.readlines()
    csc_keylist = [key[:-1] for key in csc_keys]
    with open('gmarshal_stats_feature_names.txt','r') as f:
        ptf_keys = f.readlines()
    ptf_keylist = [key[:-1] for key in ptf_keys]

    query = db.query("SELECT %s %s s.nbestobs, o.ra, o.dec FROM chandra_csc as c, ptf_obj as o, gmarshal_features_sources_R(matchedto=o,nmax=1,dmax=3) as s WHERE s.nbestobs >= 10"
        % ('c.'+'c.'.join(csc_keylist), 's.'+'s.'.join(ptf_keylist)))

    res = query.fetch()
    
    # from records is much slower than the regular df constructor
    df = pandas.DataFrame.from_records(res.as_ndarray(),index='c.csc_name')
    
    output = open('csc-ptf_features.pkl', 'wb')
    pickle.dump(df, output, pickle.HIGHEST_PROTOCOL)
    output.close()

def dump_csc_features():
    """Write CSC-only features to a pickle file"""
    db = DB(os.environ['LSD_DB'])

    with open('csc_feature_names.txt','r') as f:
        keys = f.readlines()
    keylist = [key[:-1] for key in keys]

    # this turns out not to be so slow: 20 seconds
    query = db.query("SELECT %s FROM chandra_csc" 
        % (''.join(keylist)[:-1]))

    res = query.fetch()
    
    # from records is much slower than the regular df constructor
    df = pandas.DataFrame.from_records(res.as_ndarray(),index='csc_name')
    
    output = open('csc_features.pkl', 'wb')
    pickle.dump(df, output, pickle.HIGHEST_PROTOCOL)
    output.close()

def load_csc_features(feature_filename='csc_features.pkl'):

    output = open(feature_filename, 'rb')
    df = pickle.load(output)
    output.close()

    labels = pandas.read_table('xrb_labels_csc_simbad.txt',index_col=0,
        names=['csc_name','xrb_label'])

    df = df.join(labels)

    return df

def train_rf(X, Y):

    n_folds = 5

    kfold = cross_validation.StratifiedKFold(Y, n_folds=n_folds)
    classifer = ensemble.RandomForestClassifier(n_estimators=500,class_weight='subsample',
        n_jobs=-1)
    #classifer = semi_supervised.LabelPropagation(kernel='knn')

    mean_tpr = 0.0
    mean_fpr = N.linspace(0, 1, 100)
    all_tpr = []
    precisions = N.zeros(n_folds)
    recalls = N.zeros(n_folds)
    f1s = N.zeros(n_folds)

    for i, (train, test) in enumerate(kfold):
        rf = classifer.fit(X[train],Y[train])
        probas = rf.predict_proba(X[test])
        # probas[:,1] selects probabilities that each sample is in class 1
        #fpr, tpr, thresholds = sklearn.metrics.roc_curve(Y[test], probas[:, 1])
        fpr, tpr, thresholds = sklearn.metrics.precision_recall_curve(Y[test], probas[:, 1])
        mean_tpr += scipy.interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
        roc_auc = sklearn.metrics.auc(fpr, tpr, reorder=True)
        #plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.2f)' % (i, roc_auc))

        Y_pred = rf.predict(X[test])
        precision, recall, f1, support = \
            sklearn.metrics.precision_recall_fscore_support(Y[test],Y_pred,
            beta = 1)
        print 'Fold ',i
        print 'Support: ', support
        print 'Precision: ',precision[1] # only worry about positive class
        print 'Recall: ',recall[1]
        print 'f1: ',f1[1]
        print 'Area under ROC: ', sklearn.metrics.roc_auc_score(Y[test], probas[:, 1])
        print 
        precisions[i] = precision[1]
        recalls[i] = recall[1]
        f1s[i] = f1[1]

    #plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6), label='Luck')

    mean_tpr /= len(kfold)
    mean_tpr[-1] = 1.0
    mean_auc = sklearn.metrics.auc(mean_fpr, mean_tpr, reorder=True)
    #plt.plot(mean_fpr, mean_tpr, 'k--',
    #        label='Mean ROC (area = %0.2f)' % mean_auc, lw=2)

    #plt.xlim([0.0, 0.5])
    #plt.ylim([-0.05, 1.05])
    #plt.xlabel('False Positive Rate')
    #plt.ylabel('True Positive Rate')
    #plt.title('Receiver operating characteristic')
    #plt.legend(loc="lower right")
    #plt.show()

#    importance = N.array(zip(df.columns[:-1],rf.feature_importances_),
#        dtype=[('feature','S20'),('importance','f')])
#    importance.sort(order='importance')

    print 'Mean precision: ',precisions.mean()
    print 'Mean recall: ',recalls.mean()
    print 'Mean f1: ',f1s.mean()

    return rf

def predict_rf():
    df = load_csc_features()
    labelled_df = df[df['xrb_label'] >= 0]
    unlabelled_df = df[df['xrb_label'] < 0]
    X_train = labelled_df.as_matrix()[:,:-1]
    Y_train = labelled_df.as_matrix()[:,-1]
    X_unlabelled = unlabelled_df.as_matrix()[:,:-1]

    classifer = ensemble.RandomForestClassifier(n_estimators=500,
        n_jobs=16,compute_importances=True)

    # learn on all labelled features
    rf = classifer.fit(X_train,Y_train)
    probas = rf.predict_proba(X_unlabelled)

    # probas[:,1] selects probabilities that each sample is in class 1
    probas_series = pandas.Series(probas[:,1],index=unlabelled_df.index,
        name='probability')

    result_df = unlabelled_df.join(probas_series)
    str = result_df.to_string(columns=['ra','dec','probability'],col_space=1)
    with open('binary_probability_unlabelled_csc.txt','w') as f:
        f.write(str)
        f.close()

