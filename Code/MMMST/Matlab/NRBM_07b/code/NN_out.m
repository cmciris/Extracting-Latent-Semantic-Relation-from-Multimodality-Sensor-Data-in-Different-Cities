function Xout = NN_out(nn,data,K)
%function Xout = NN_out(nn,data,K)
% Feed the data to the neural net and get the output at the layer k
% X : [d n]
% nn : trained neural net
% e.g. nn =
%      numhids: [200 200 2]
%    hidbiases: {3x1 cell}
%       vishid: {3x1 cell}

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

if isempty(data)
	Xout = [];
	return;
end

numhids = nn.numhids;
hidbiases = nn.hidbiases;
vishid = nn.vishid;

[be en] = make_minibatches(size(data,2),1000);

if (nargin<3) K = numel(numhids);end

DIM = size(vishid{K},2);
Xout = zeros(DIM,size(data,2));
for b=1:numel(be)
	X = data(:,be(b):en(b));
	n = size(X,2);
	for k=1:K
	  X = 1./(1 + exp(-vishid{k}'*X - repmat(hidbiases{k}',1,n)));
	end
	Xout(:,be(b):en(b)) = X;
end
