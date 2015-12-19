
date


####################################################################################################################################
# 
# Executing the Random Forest Classifier to train for RR Lyrae and
# QSO candidate classification
# 
# 
####################################################################################################################################


PYTHONCOMPILED=MODULE
export PYTHONCOMPILED

export NWORKERS=16
export OMP_NUM_THREADS=1


echo $!

#train
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 1 &#
echo $!
let PID=$!
wait $PID

date  

echo "training RR Lyrae done"


python RF_classification_DRW_grid_QSO_alldata_noEBV_noSDSS.py 1 &

echo $!
let PID=$!
wait $PID

echo "training QSO done"




date  

 
