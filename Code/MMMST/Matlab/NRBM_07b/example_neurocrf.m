if ~exist('PATH_OK')
  addpath code; % path to NRBM and CRF libraries
  addpath data; % path to ocr data
  addpath(genpath('soft/stprtool'));
  PATH_OK = 1;
end

SMALL_SETTING = 1;

fold = 1; % select fold number in cross-validation
arch = [200 200];lrate = 1e-1; maxepoch = 50; % neural net parameters
criterion = 'CML'; regmod = 'uniform'; % neurocrf learning parameters

if SMALL_SETTING
	trainsize = 1; lambda = 1e-4; % this is a good value for small setting
else
	trainsize = 9; lambda = 2e-5; % this is a good value for large setting
end

[trainData testData] = getdata_ocr('data/ocr.mat',fold,trainsize);

neurocrf0	= NNstruct_init(arch,trainData,'lambda',lambda,'lrate',lrate,'maxepoch',maxepoch);
neurocrf	= NNstruct_train(neurocrf0,trainData,'lambda',lambda,'flag',4,'criterion',criterion,'regmod',regmod);

ErrTrain = NNstruct_eval(neurocrf,trainData);
ErrTest = NNstruct_eval(neurocrf,testData);
disp('===============================================================');
disp(' Neural conditional random field');
disp('===============================================================');
disp(sprintf(' Neural net architecture : [ %s]',sprintf('%d ',arch)));
disp(sprintf(' Discriminative training criterion : %s',criterion));
disp(sprintf(' Regularization parameter (lambda) : %.0e',lambda));
disp(sprintf(' Cross-validation-fold=%d trainsize=%d fold(s)', fold,trainsize));
disp(sprintf(' ErrTrain=%.4f   ErrTest=%.4f',ErrTrain,ErrTest));
disp('===============================================================');

% save my_neurocrf.mat neurocrf  (?)
