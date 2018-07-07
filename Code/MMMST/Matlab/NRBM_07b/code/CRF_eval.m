function error_rate = CRF_eval(model,data)
%function error_rate = CRF_eval(model,data)
%
% work with dense data
% e.g. trainData =
%		       X: [128x4617 double]
%		       Y: [4617x1 double]
%		     Seq: [626x2 double]

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

	numLab = model.numLab;
	w = model.w;
	dim = (numel(w)-numLab*numLab)/numLab;

	lab = mexViterbiInference(data.X,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);

	error_rate = sum(lab(:)~=data.Y) / numel(lab);

