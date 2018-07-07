cd soft/stprtool
compilemex

cd ../../code 
mex wsum_col.c 
mex wsum_row.c
mex mexForwardBackward.c 
mex mexForwardBackwardBp.c
mex mexPosLearnMostViolationBp.c
mex mexViterbiMostViolation.c
mex mexViterbiInference.c
mex mexViterbiMostViolationBp.c
mex mexPosLearnMostViolation.c
mex mexViterbiInferenceSegmented.c

cd ../
