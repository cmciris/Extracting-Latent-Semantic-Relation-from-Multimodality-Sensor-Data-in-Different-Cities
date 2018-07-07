function neurocrf = NNstruct_train(neurocrf0,trainData,varargin)
%function neurocrf = NNstruct_train(neurocrf0,trainData,varargin)
%
% option : lambda(1.0), maxiter (500), maxCP (200), epsilon(0.01), flag(0) , testData([])
%
% work with dense data
% e.g. trainData =
%		       X: [4617x128 double]
%		       Y: [4617x1 double]
%		     Seq: [626x2 double]
%

%
% Copyright (C) 2012, by Trinh-Minh-Tri Do minhtrido@gmail.com
%
%   This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get parameter
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	lambda  	= get_varargin( varargin,  'lambda',    1);	% regularization parameter
	maxiter 	= get_varargin( varargin, 'maxiter',  500);
	maxCP   	= get_varargin( varargin,   'maxCP',  200);	% maximum number of cutting plane
	epsilon 	= get_varargin( varargin, 'epsilon', 0.01);
	flag		= get_varargin( varargin,    'flag', 0);
	testData	= get_varargin( varargin,'testData', []);
	criterion	= get_varargin( varargin,'criterion', 'CML'); % CML : conditional likelihood, 
	regmod		= get_varargin( varargin,'regmod','uniform');   % recommendation : uniform or stdcrf
	gsd			= get_varargin( varargin,'gsd',1);					% regularization coefficient for 1st layer
	fullfeature	= get_varargin( varargin,'fullfeature',0);			% say if we use feature from all layers or only last layer
	lambdaTrans = get_varargin( varargin, 'lambdaTrans',            1);	% regularization coefficient for transition probability
	bigram		= get_varargin( varargin,      'bigram',           []);	% bigram to be regularized
	neurocrfWarmStart	= get_varargin( varargin,    'neurocrfWarmStart',  []);	% warm start solution
	minibatchsize		= get_varargin( varargin,    'minibatchsize',     100);	% the size minibatch


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Initial solution
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(neurocrfWarmStart)
		w0	= [NN_getw(neurocrf0.nn);neurocrf0.crf.w(:)];
	else
		w0	= [NN_getw(neurocrfWarmStart.nn);neurocrfWarmStart.crf.w(:)];
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Prepair
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	numLab	= max(trainData.Y);
	dim     = size(trainData.X,1);
	nn      = neurocrf0.nn;
	crf	= neurocrf0.crf;
	dim_crf = numel(crf.w);
        dim_nn  = numel(w0)-dim_crf;
	dim_crf_trans	= numLab*numLab;
	dim_crf_local	= dim_crf - dim_crf_trans;
	dim_nn_VH1 = numel(nn.vishid{1}) + numel(nn.hidbiases{1});

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Initial solution
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if isempty(bigram)
		bigram = zeros(numLab,numLab);
	end
	wreg	= [NN_getw(nn);bigram(:);zeros(dim_crf_local,1)];
	if strcmp(regmod,'uniform')
        reg     = [1*ones(dim_nn_VH1,1); 1*ones(dim_nn-dim_nn_VH1,1); lambdaTrans*ones(dim_crf_trans,1); 1*ones(dim_crf_local,1)];
	elseif strcmp(regmod,'normalized')
		mCRF	= mean(abs(crf.w));	% mean value of last layer CRF
		mCRF	= max(mCRF,1e-2);
		K	= numel(nn.vishid);	% for each layer
		reg	= [];
		for k=1:K
			mNN = mean(abs(nn.vishid{k}(:)));
			reg = [reg; (mCRF/mNN)*ones(numel(nn.vishid{k}),1)];
			mNN = mean(abs(nn.hidbiases{k}(:)));
			reg = [reg; (mCRF/mNN)*ones(numel(nn.hidbiases{k}),1)];
		end
		reg = [reg;lambdaTrans*ones(dim_crf_trans,1);1*ones(dim_crf_local,1)];
	elseif strcmp(regmod,'stdone')
		MIN_STD = 0.01;
		K	= numel(nn.vishid);	% for each layer
		reg	= [];
		for k=1:K
			stdNN = max(std(nn.vishid{k}(:)),MIN_STD);
			reg = [reg;(1/stdNN)*ones(numel(nn.vishid{k}),1)];
			stdNN = max(std(nn.hidbiases{k}(:)),MIN_STD);
			reg = [reg;(1/stdNN)*ones(numel(nn.hidbiases{k}),1)];
		end
		reg = [reg;lambdaTrans*ones(1,dim_crf_trans);1*ones(dim_crf_local,1)];
	elseif strcmp(regmod,'stdcrf')
		MIN_STD = 0.01;
		stdcrf  = max(std(crf.w((dim_crf_trans+1):end)),MIN_STD);
		K	= numel(nn.vishid);	% for each layer
		reg	= [];
		for k=1:K
			stdNN = max(std(nn.vishid{k}(:)),MIN_STD);
			reg = [reg;(stdcrf/stdNN)*ones(numel(nn.vishid{k}),1)];
			stdNN = max(std(nn.hidbiases{k}(:)),MIN_STD);
			reg = [reg;(stdcrf/stdNN)*ones(numel(nn.hidbiases{k}),1)];
		end
		reg = [reg;lambdaTrans*ones(dim_crf_trans,1); 1*ones(dim_crf_local,1)];
	elseif strcmp(regmod,'gsd')
        	reg     = [gsd*ones(dim_nn_VH1,1);1*ones(dim_nn-dim_nn_VH1,1);lambdaTrans*ones(dim_crf_trans,1);1*ones(dim_crf_local,1)];
	else
		disp(sprintf('\tUnknown regmod option ''%s''. Allowed value : uniform,normalized,stdone,stdcrf,gsd',regmod));
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Preparing call back function (risk and its subgradient) and its required data
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	neurocrf = neurocrf0; 				% copying structure of neurocrf0
	trainDataCell = split_data(trainData,'minibatchsize',minibatchsize);	% Avoid using too much memory by spliting training data in minibatch of sequences.
								% This can be used for parallel computing also.
								% If parfor is supported then just modify the for loop in FGrad.
	auxdata = {trainData trainDataCell testData dim neurocrf flag criterion fullfeature};
	handleFGrad	= @FGrad;
	handleFReport = @test_model;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calling optimization routine : Generalized Non-convex Regurarized Bundle Method
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	[w stats] = NRBM(handleFGrad,auxdata,lambda,w0,...
				'wreg',wreg,...
				'reg',reg,...
				'freport',handleFReport,...
				'computeGapQP',0,'verbosity',flag>0,'Rconvex',0,...
				'maxiter',maxiter,...
				'maxCP',maxCP,...
				'epsilon',epsilon,'LS',1);

	neurocrf.nn = NN_setw(neurocrf.nn,w(1:dim_nn));
	neurocrf.crf.w = w((dim_nn+1):end);

%=======================================================
function [F,Grad,nothing] = FGrad(w,auxdata)
	
	[data datacell testData dim neurocrf flag criterion fullfeature] = deal(auxdata{:});

	%%% set new parameters to neurocrf %%%
	neurocrf.nn  = NN_setw(neurocrf.nn,w);
	dim_nn       = NN_getdim(neurocrf.nn);
	neurocrf.crf.w = w((dim_nn+1):end);

	F = 0;
	Grad = zeros(numel(w),1);
	n = 0;
	for i=1:numel(datacell)
		if fullfeature
		        [f grad_neurocrf] = FGradNeuroStructFullFeature(neurocrf,datacell{i},criterion);
		else
		        [f grad_neurocrf] = FGradNeuroStruct(neurocrf,datacell{i},criterion);
		end
		F	= F + f;
		Grad	= Grad + [NN_getw(grad_neurocrf.nn);grad_neurocrf.crf.w];
		n	= n + size(datacell{i}.X,2);
	end
	
	
	F = F/n;
	Grad = Grad / n;
	nothing = [];

%	keyboard;

%=======================================================
function str=test_model(w,auxdata)
	[trainData datacell testData dim neurocrf flag criterion fullfeature] = deal(auxdata{:});
	neurocrf.nn  = NN_setw(neurocrf.nn,w);
	dim_nn       = NN_getdim(neurocrf.nn);
	neurocrf.crf.w = w((dim_nn+1):end);
	str = ' ';
	if flag==1
		err_train = NNstruct_eval(neurocrf,trainData,'fullfeature',fullfeature);
		str = sprintf('ErrTrain=%.4f ',err_train);
	end
	if (flag==2)
		if isempty(testData)
			err_test = 1;
		else
			err_test = NNstruct_eval(neurocrf,testData,'fullfeature',fullfeature);
		end
		str = sprintf('ErrTest=%.4f ',err_test);
	end
	if flag==3
		err_train = NNstruct_eval(neurocrf,trainData,'fullfeature',fullfeature);
		if isempty(testData)
			err_test = 1;
		else
			err_test = NNstruct_eval(neurocrf,testData);
		end
		str = sprintf('ErrTrain=%.4f ErrTest=%.4f',err_train,err_test);
	end

%===========================================================
function [f grad_neurocrf] = FGradNeuroStruct(neurocrf,data,criterion)
  nn = neurocrf.nn;
  crf = neurocrf.crf;
  K = numel(nn.numhids);
  grad_neurocrf = neurocrf; % just want to copy the structure

  % feed forward neural
  X = NN_forward(nn,data.X); % [v n]

  % compute marginal probability
  datacrf = data;
  datacrf.X = X{K+1};
  [f grad_crf bp] = Marginal(crf,datacrf,criterion); % bp : [dim(k+1) n] matrix

  grad_neurocrf.crf.w = grad_crf(:);
  nbtok = size(data.X,2);
  ONE = ones(1,nbtok);

  % back propagation
  for k=K:-1:1
    w = [nn.vishid{k};nn.hidbiases{k}];% [ dim(k) x dim(k+1) ]
    bp = bp .* X{k+1} .* (1-X{k+1});   % [ dim(k+1) n]
    dw = [X{k};ONE]*bp';               % [ dim(k) x     ntok ] * [ ntok     x dim(k+1)]
    bp = w * bp ;                      %  [ dim(k) x dim(k+1) ] * [ dim(k+1) n]
    bp = bp(1:end-1,:);
    grad_neurocrf.nn.vishid{k}    = dw(1:end-1,:);
    grad_neurocrf.nn.hidbiases{k} = dw(end,:);
  end

%===========================================================
function [f grad_neurocrf] = FGradNeuroStructFullFeature(neurocrf,data,criterion)
  nn = neurocrf.nn;
  crf = neurocrf.crf;
  K = numel(nn.numhids);
  grad_neurocrf = neurocrf; % just want to copy the structure

  % feed forward neural
  X = NN_forward(nn,data.X); % [v n]

  % compute marginal probability
  datacrf = data;
  for k=1:K
    datacrf.X = [datacrf.X;X{k+1}];
  end
  [f grad_crf BP] = Marginal(crf,datacrf,criterion); % BP : [(dim(V)+dim(H1)+dim(H2)+...+dim(HK)) n] matrix

  grad_neurocrf.crf.w = grad_crf;
  nbtok = size(data.X,2);
  ONE = ones(1,nbtok);
  v = size(data.X,1);

  en = cumsum([v nn.numhids]);be = [1 (en(1:end-1)+1)];

  % back propagation
  bp = BP(be(K+1):en(K+1),:);
  for k=K:-1:1
    w = [nn.vishid{k};nn.hidbiases{k}];% [ dim(k) x dim(k+1) ]
    bp = bp .* X{k+1} .* (1-X{k+1});   % [ dim(k+1) n]
    dw = [X{k};ONE]*bp';               % [ dim(k) x     ntok ] * [ ntok     x dim(k+1)]
    bp = w * bp ;                      %  [ dim(k) x dim(k+1) ] * [ dim(k+1) n]
    bp = bp(1:end-1,:) + BP(be(k):en(k),:);
    grad_neurocrf.nn.vishid{k}    = dw(1:end-1,:);
    grad_neurocrf.nn.hidbiases{k} = dw(end,:);
  end

%===========================================================
function [f grad_crf bp] = Marginal(crf,data,criterion)
  if strcmp(criterion,'CML')
	  [f grad_crf bp] = mexForwardBackwardBp(...
			data.X,...
			data.Y,...
			data.Seq(:,1),data.Seq(:,2),...
			crf.numLab,crf.dim,...
			crf.w);	
  elseif strcmp(criterion,'M3N')
	  [f grad_crf bp] = mexViterbiMostViolationBp(...
			data.X,...
			data.Y,...
			data.Seq(:,1),data.Seq(:,2),...
			crf.numLab,crf.dim,...
			crf.w);	
  elseif strcmp(criterion,'PosLearn')
	  [f grad_crf bp] = mexPosLearnMostViolationBp(...
			data.X,...
			data.Y,...
			data.Seq(:,1),data.Seq(:,2),...
			crf.numLab,crf.dim,...
			crf.w);
  else
	disp(sprintf('NNstruct_train.m : unknown learning criterion %s',criterion));update_diary();
  end

%===========================================================
function datacell = split_data(data,varargin)
% e.g. data =
%		       X: [128x4617 double]
%		       Y: [4617x1 double]
%		     Seq: [626x2 double]

minibatchsize	= get_varargin( varargin,'minibatchsize',100);

n		= size(data.Seq,1);
numbatches	= -floor(-n/minibatchsize);
be 		= round((0:(numbatches-1))*(n/numbatches))+1;
en 		= [be(2:end)-1,n];

datacell = cell(numbatches,1);
for i=1:numbatches
	be1	= data.Seq(be(i),1);
	en1	= data.Seq(en(i),2);
	datacell{i}.X = data.X(:,be1:en1);
	datacell{i}.Y = data.Y(be1:en1);
	datacell{i}.Seq = data.Seq(be(i):en(i),:) - (be1-1);
end
