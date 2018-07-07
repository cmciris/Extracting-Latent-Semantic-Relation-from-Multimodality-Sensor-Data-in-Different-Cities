function dim_nn = NN_getdim(nn)
% function dim_nn = NN_getdim(nn)
%
% get number of parameters of neural network
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

  numhids = nn.numhids;
  hidbiases = nn.hidbiases;
  vishid = nn.vishid;

  K = numel(numhids);

  dim_nn = 0;
  for k=1:K
    dim_nn = dim_nn + numel(nn.vishid{k}) + numel(nn.hidbiases{k});
  end
