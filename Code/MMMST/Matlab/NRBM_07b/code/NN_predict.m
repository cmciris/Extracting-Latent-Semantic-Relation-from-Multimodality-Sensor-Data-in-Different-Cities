function targetout = NN_predict(nn,data)
%function targetout = NN_predict(nn,data)
% nn	  : neural net
% data    : [NxD]
% output  : [NxM]

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

	w1=[nn.vishid{1}; nn.hidbiases{1}];
	w2=[nn.vishid{2}; nn.hidbiases{2}];
	w3=[nn.vishid{3}; nn.hidbiases{3}];
	w_class = nn.w_class;
    N = size(data,1);
    M = size(w_class,2);
	l1=size(w1,1)-1;
	l2=size(w2,1)-1;
	l3=size(w3,1)-1;
	l4=size(w_class,1)-1;
	l5=M; 

	XX = [data ones(N,1)];
	w1probs = 1./(1 + exp(-XX*w1)); w1probs = [w1probs  ones(N,1)];
	w2probs = 1./(1 + exp(-w1probs*w2)); w2probs = [w2probs ones(N,1)];
	w3probs = 1./(1 + exp(-w2probs*w3)); w3probs = [w3probs  ones(N,1)];
	targetout = 1./(1 + exp(-w3probs*w_class));

