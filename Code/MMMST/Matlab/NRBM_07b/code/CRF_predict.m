function Ypredict = CRF_predict(model,data,varargin)
%function Ypredict = CRF_predict(model,data,varargin)
%
% work with dense data
% e.g. trainData =
%		       X: [128x4617 double]
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

	Q = get_varargin(varargin, 'Q',[]);

	numLab = model.numLab;
	w = model.w;
	dim = (numel(w)-numLab*numLab)/numLab;

	if isempty(Q)
		Ypredict = mexViterbiInference(data.X,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);
	else
		Ypredict = mexViterbiInferenceSegmented(data.X,Q,...
					data.Seq(:,1),data.Seq(:,2),...
					numLab,dim,...
					w);
	end
