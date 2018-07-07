if ~exist('PATH_OK')
  addpath code; % path to NRBM and CRF libraries
  addpath data; % path to ocr data
  addpath(genpath('soft/stprtool'));
  PATH_OK = 1;
end

lambdaCML = 5e-4;
lambdaM3N = 1e-3;
lambdaPosLearn = 1e-2;

disp(' training CRF with M3N ...'); 
[errTrainM3N errTestM3N] = run_crf_ocr('fold',1,'trainsize',1,'lambda',lambdaM3N,'criterion','M3N','epsilon',0.001);

if 0
disp(' training CRF with CML ...'); 
[errTrainCML errTestCML] = run_crf_ocr('fold',1,'trainsize',1,'lambda',lambdaCML,'criterion','CML');
disp(' training CRF with M3N ...'); 
[errTrainM3N errTestM3N] = run_crf_ocr('fold',1,'trainsize',1,'lambda',lambdaM3N,'criterion','M3N');
disp(' training CRF with PosLearn ...'); 
[errTrainPosLearn errTestPosLearn] = run_crf_ocr('fold',1,'trainsize',1,'lambda',lambdaPosLearn,'criterion','PosLearn');

disp(''); 
disp('============================================================================================='); 
disp('=                      Final results of CRFs with different criteria                        '); 
disp('============================================================================================='); 
disp(sprintf('= Conditional Maximum Likelihood (CML) lambda = %.1e : errTrain=%.4f errTest=%.4f',lambdaCML,errTrainCML,errTestCML)); 
disp(sprintf('= Max-margin                     (M3N) lambda = %.1e : errTrain=%.4f errTest=%.4f',lambdaM3N,errTrainM3N,errTestM3N)); 
disp(sprintf('= Max-margin Position Learn (PosLearn) lambda = %.1e : errTrain=%.4f errTest=%.4f',lambdaPosLearn,errTrainPosLearn,errTestPosLearn)); 
disp('============================================================================================='); 
end
