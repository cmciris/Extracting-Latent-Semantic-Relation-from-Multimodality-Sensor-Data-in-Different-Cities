function [ErrTrain ErrTest]= run_crf_ocr(varargin)
%function [ErrTrain ErrTest] = run_crf_ocr(varargin)
%
% Example [ErrTrain ErrTest] = run_crf_ocr('fold',1,'trainsize',9,'lambda',1e-4);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% debug section 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	FORCE_LEARN_NEUROCRF	=	0;

	fdata		= get_varargin(varargin,     'fdata', 'ocr.mat');
	fold		= get_varargin(varargin,      'fold',         1);
	trainsize	= get_varargin(varargin, 'trainsize',         1);
	lambda		= get_varargin(varargin,    'lambda',      1e-4);
	criterion	= get_varargin(varargin, 'criterion',      'CML');
	epsilon		= get_varargin(varargin,   'epsilon',       1e-2);
	flag	= get_varargin(varargin,    'flag',    1);
	
	[trainData testData] = getdata_ocr(fdata,fold,trainsize);
	
%  	crf = CRF_train(trainData,'lambda',lambda,'flag',flag,'criterion',criterion,'epsilon',epsilon,'testData',testData); % with test data
	crf = CRF_train(trainData,'lambda',lambda,'flag',flag,'criterion',criterion,'epsilon',epsilon); % no test
	ErrTrain = CRF_eval(crf,trainData);
	ErrTest = CRF_eval(crf,testData);
	disp(sprintf('*** CRF-%s  linear %.0e fold%2d fold in train %2d *** ErrTrain  %.4f   ErrTest  %.4f', criterion,lambda,fold,trainsize,ErrTrain,ErrTest));

