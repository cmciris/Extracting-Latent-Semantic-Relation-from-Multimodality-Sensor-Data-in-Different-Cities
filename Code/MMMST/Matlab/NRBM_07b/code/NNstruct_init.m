function neurocrf = NNstruct_init(arch,trainData,varargin)
%function neurocrf = NNstruct_init(arch,trainData,varargin)
%
% Initialize NNstruct with RBM then CRF
%
% arch		: number of hidden unit in each layer
% trainData	: fully labelled sequential data , work with dense data
% 		  e.g. trainData =
%		       X: [128x4617 double]
%		       Y: [4617x1 double]
%		     Seq: [626x2 double]
%
% additional option : lambda(1e-4), nn0([]) , maxepoch(50), batchsize(100)
%
% Example : NNstruct_init([200 200 50],trainData,'lambda',2e-4,'nn',nn);
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

	nn			= get_varargin(varargin,            'nn',    []); % a separately learned neural net can be used
	verbosity	= get_varargin(varargin,     'verbosity',     1);

	% Neural network part
	maxepoch	= get_varargin(varargin,      'maxepoch',    50);
	batchsize	= get_varargin(varargin,     'batchsize',   100);
	lrate		= get_varargin(varargin,      'lrate',           1e-1);   % lrate
	gaussian	= get_varargin(varargin,      'gaussian',     0);

	% CRFs part
	lambda  	= get_varargin(varargin,        'lambda',  1e-4);
	criterion	= get_varargin(varargin,     'criterion',  'CML');
	fullfeature	= get_varargin(varargin,   'fullfeature',     0);
	lambdaTrans	= get_varargin(varargin, 'lambdaTrans',            1);   % regularization coefficient for transition probability
	bigram		= get_varargin(varargin,      'bigram',           []);   % a bigram can be used for regularization here


	%%%%%%%%%%%%%%%%%%%%%%%%%
	% init Neural network
	%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty(nn)
		nn = NN_init(trainData.X,arch,...
				'lrate',lrate,...
				'maxepoch',maxepoch,...
				'batchsize',batchsize,...
				'verbosity',verbosity,...
				'gaussian',gaussian);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% init last layer with structured model CRF or M3N
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	trainDataT = trainData;
	if fullfeature
		K = numel(nn.numhids);
		X = NN_forward(nn,trainData.X); % [v n]
		for k=1:K
			trainDataT.X = [trainDataT.X;X{k+1}];
		end    
		clear X;
	else
		trainDataT.X = NN_out(nn,trainData.X);
	end
	crf = CRF_train(trainDataT,'lambda',lambda,'flag',4,'lambdaTrans',lambdaTrans,'bigram',bigram,'criterion',criterion); % notest

	neurocrf = struct('nn',nn,'crf',crf);
