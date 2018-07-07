function Ypredict = NNstruct_predict(neurocrf,data,varargin)
%function Ypredict = NNstruct_predict(neurocrf,data,varargin)

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

  fullfeature	= get_varargin(varargin, 'fullfeature', 0);   

  nn = neurocrf.nn;
  crf = neurocrf.crf;

  datacell = split_data(data);
  Ypredict = [];
  for k=1:numel(datacell)
	datacrf = datacell{k};
	if fullfeature
		datacrf.X = NN_outFF(nn,datacrf.X);
	else
		datacrf.X = NN_out(nn,datacrf.X);
	end
	Ypredict = [Ypredict;M3N_predict(crf,datacrf)];
  end

%===========================================================
function datacell = split_data(data,varargin)
% e.g. data =
%		       X: [128x4617 double]
%		       Y: [4617x1 double]
%		     Seq: [626x2 double]

batchsize	= get_varargin(varargin,'batchsize',100);

n		= size(data.Seq,1);
numbatches	= -floor(-n/batchsize);
be 		= round((0:(numbatches-1))*(n/numbatches))+1;
en 		= [be(2:end)-1,n];

datacell = cell(numbatches,1);
for i=1:numbatches
	be1	= data.Seq(be(i),1);
	en1	= data.Seq(en(i),2);
	datacell{i}.X = data.X(:,be1:en1);
	datacell{i}.Y = data.Y(be1:en1);
	datacell{i}.Seq = data.Seq(be(i):en(i),:) - (be1-1);
end
