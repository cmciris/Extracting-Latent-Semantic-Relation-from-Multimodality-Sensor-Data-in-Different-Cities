function X = NN_forward(nn,data)
%function X = NN_forward(nn,data)
% Forward neural network for data and get output at each layer (input included)
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

  numhids   = nn.numhids;
  hidbiases = nn.hidbiases;
  vishid    = nn.vishid;
  K         = numel(numhids);

  X         = cell(K+1,1);
  X{1}      = data;
  n         = size(data,2);
  for k=1:K
    X{k+1}  = 1./(1 + exp(-vishid{k}'*X{k} - repmat(hidbiases{k}',1,n)));
  end
