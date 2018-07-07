function nn_out = NN_flip(nn)
%function nn_out = NN_flip(nn)

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

nn_out = struct('numhids',nn.numhids);
K = numel(nn.numhids);
nn_out.hidbiases = cell(K,1);
nn_out.vishid = cell(K,1);
nn_out.visbiases = cell(K,1); 

for k =1:K
	nn_out.hidbiases(k) = nn.visbiases(K-k+1);
	nn_out.vishid{k} = nn.vishid{k}';
	nn_out.visbiases(k) = nn.hidbiases(K-k+1);
end
