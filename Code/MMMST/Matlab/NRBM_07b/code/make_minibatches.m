function [be en] = make_minibatches(n,batchsize)
%function [be en] = make_minibatches(n,batchsize)
%
% n           number of element
% batchsize   size of minibatche

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

numbatches	= -floor(-n/batchsize);
be 		= round((0:(numbatches-1))*(n/numbatches))+1;
en 		= [be(2:end)-1,n];
