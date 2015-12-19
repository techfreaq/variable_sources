
date

####################################################################################################################################
# 
# Random Forest Classification for RR Lyrae candidates, applied for the patches shown below
# 
# 
# 
####################################################################################################################################

PYTHONCOMPILED=MODULE
export PYTHONCOMPILED

export NWORKERS=16
export OMP_NUM_THREADS=1


echo $!



python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch1/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch1/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID

date  

echo "patch 1 done"




python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch10/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch10/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID

date  

echo "patch 10 done"


 
 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch2/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch2/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &
 
echo $!
let PID=$!
wait $PID
 
date  
 
echo "patch 2 done"
 
 
 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch3/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch3/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID
 
date  
 
echo "patch 3 done"
 
 
 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch4/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch4/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &
 
echo $!
let PID=$!
wait $PID
 
date  
 
echo "patch 4 done"
 
 
 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch5/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch5/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID

date  

echo "patch 5 done"
 
 
 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch6/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch6/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &
 
echo $!
let PID=$!
wait $PID

date  

echo "patch 6 done"



python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch7/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch7/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID

date  

echo "patch 7 done"


python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch9/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch9/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID

date  

echo "patch 9 done"
 
 


python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch11/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch11/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &

echo $!
let PID=$!
wait $PID
 
date  
 
echo "patch 11 done"
 

 
python RF_classification_DRW_grid_RRLyrae_alldata_noEBV_noSDSS.py 0 "/a41233d1/hernitschek/_results_juelich/_results_strucfunc_new/_patch12/" "/a41233d1/hernitschek/_results_juelich/_results_classifier/results_pRRLyrae_new_noSDSS/_patch12/" "/a41233d1/hernitschek/_sdssdata_targetsets_newsdss/" &
 
echo $!
let PID=$!
wait $PID

date  

echo "patch 12 done"



date  

 
