function neurocrf_scale = NNstruct_scaleTrans(neurocrf,scaleDiag,scaleJump);
%function neurocrf_scale = NNstruct_scaleTrans(neurocrf,scaleDiag,scaleJump);

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

    numLab	= neurocrf.crf.numLab;
	dimTrans = numLab*numLab;
	matscale = ones(numLab,numLab)*scaleJump;
	for y=1:numLab
		matscale(y,y) = scaleDiag;
	end
	neurocrf_scale = neurocrf;
	neurocrf_scale.crf.w(1:dimTrans) = neurocrf.crf.w(1:dimTrans) .* reshape(matscale,1,dimTrans);
