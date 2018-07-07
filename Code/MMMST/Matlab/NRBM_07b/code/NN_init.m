function nn = NN_init(X,numhids,varargin)
% function nn = NN_init(X,numhids,varargin)
% init neural network with rbm
% X : [DxN]
% numhids : [1xK] where K is the number of layers, and numhids(k) is number of hidden nodes in layer k
% 
% options : maxepoch(50) batchsize(100) verbosity(1) gaussian(0)
% Example : nn = NN_init(randn(1000,20),[30 50 30],'gaussian',1);
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

maxepoch	= get_varargin(varargin,  'maxepoch',  100);
batchsize	= get_varargin(varargin, 'batchsize',  100);
verbosity	= get_varargin(varargin, 'verbosity',    1);
gaussian	= get_varargin(varargin, 'gaussian',     0);
lrate		= get_varargin(varargin, 'lrate',     1e-1);
rbm_init	= get_varargin(varargin,'rbm_init','randn');
D		= get_varargin(varargin,'D',3);
momentum	= get_varargin(varargin,'momentum',0);
gsd		= get_varargin(varargin,'gsd',1);  % set standard deviation of visible gaussian unit
XHoldOut	= get_varargin(varargin,'XHoldOut',[]); % Hold-out data [dim x num_ex]

K = numel(numhids);
vishid = cell(K,1);
hidbiases = cell(K,1);
[dim n] = size(X);
for k=1:K
	dispV(1,verbosity,sprintf('Init layer %d/%d with RBM : %d x %d',k,K,size(X,1),numhids(k)));
	[vishid{k},hidbiases{k}] = rbm(X,numhids(k),'maxepoch',maxepoch,...
		'batchsize',batchsize,'verbosity',verbosity,...
		'gaussian',(gaussian)&&(k==1),'lrate',lrate,...
		'rbm_init',rbm_init,'D',D,'momentum',momentum,...
		'gsd',gsd,'XHoldOut',XHoldOut);
	if (k<K)
		nn_this.numhids = numhids(k);
		nn_this.vishid  = vishid(k);
		nn_this.hidbiases = hidbiases(k);
		X = NN_out(nn_this,X);
		XHoldOut = NN_out(nn_this,XHoldOut);
	end
end

nn = struct('numhids',numhids);
nn.hidbiases = hidbiases;
nn.vishid = vishid;

