function model = CRF_train(trainData,varargin)
%function model = CRF_train(trainData,varargin)
%
% option : lambda(1.0), maxiter (500), maxCP (200), epsilon(0.01)
%
% work with dense data
% e.g. trainData =
%       X: [128x4617 double]
%       Y: [4617x1 double]
%     Seq: [626x2 double]

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

	lambda  = get_varargin(varargin,  'lambda', 1e-2);
	maxiter = get_varargin(varargin, 'maxiter',  500);
	maxCP   = get_varargin(varargin,   'maxCP',  200);
	epsilon = get_varargin(varargin, 'epsilon', 0.01);
	flag	= get_varargin(varargin,    'flag',    1);
	testData	= get_varargin(varargin,'testData',   []);
	criterion	= get_varargin(varargin,'criterion', 'CML');
	lambdaTrans	= get_varargin(varargin, 'lambdaTrans',            1);   % regularization coefficient for transition probability
	bigram		= get_varargin(varargin,      'bigram',           []);   % bigram to be regularized

	numLab	= max(trainData.Y);
	dim     = size(trainData.X,1);


	if isempty(bigram)
		bigram = zeros(numLab,numLab);
	end

	w0	= [bigram(:); rand(dim*numLab,1) * 0.1 * 0];
	wreg	= [bigram(:); zeros(dim*numLab,1)];
	reg	= [ones(numLab*numLab,1)*lambdaTrans;ones(dim*numLab,1)];

	handleFGrad	= @FGrad;
	handleFReport	= @test_model;


	model.w = w0;
	model.dim = dim;
	model.numLab = numLab;

	auxdata = {trainData testData dim model flag criterion};

	[w stats] = NRBM(handleFGrad,auxdata,lambda,w0,...
				'wreg',wreg,...
				'reg',reg,...
				'freport',handleFReport,...
				'computeGapQP',0,'verbosity',flag>0,'Rconvex',1,...
				'maxiter',maxiter,...
				'maxCP',maxCP,...
				'epsilon',epsilon,'LS',1);
	model.w = w;

%===========================================================================================
function [F,Grad,nothing] = FGrad(w,auxdata)
	[data testData dim model flag criterion] = deal(auxdata{:});
	numLab = model.numLab;
	if strcmp(criterion,'CML')
		[F,Grad] = mexForwardBackward(data.X,...
					data.Y,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);
	elseif strcmp(criterion,'M3N')
		[F,Grad] = mexViterbiMostViolation(data.X,...
					data.Y,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);
	elseif strcmp(criterion,'PosLearn')
		[F,Grad] = mexPosLearnMostViolation(data.X,...
					data.Y,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);
	else
		disp(sprintf('invalide criteria (%s) for learning structured model ! Allowed : CML M3N PosLearn',criteria));
		update_diary();
		return;
	end
	n = size(data.Y,1);
	F = F/n;
	Grad = Grad(:) / n;
	nothing = [];

%===========================================================================================
function str=test_model(w,auxdata)
	[trainData testData dim model flag criteria] = deal(auxdata{:});
	model.w = w;
	str = ' ';
	if flag==1
		err_train = CRF_eval(model,trainData);
		str = sprintf('ErrTrain=%.4f',err_train);
	end
	if flag==2
		if isempty(testData)
			err_test=1;
		else
			err_test = CRF_eval(model,testData);
		end
		str = sprintf('ErrTest=%.4f',err_test);
	end
	if flag==3
		err_train = CRF_eval(model,trainData);
		if isempty(testData)
			err_test=1;
		else
			err_test = CRF_eval(model,testData);
		end
		str = sprintf('ErrTrain=%.4f ErrTest=%.4f',err_train,err_test);
	end
