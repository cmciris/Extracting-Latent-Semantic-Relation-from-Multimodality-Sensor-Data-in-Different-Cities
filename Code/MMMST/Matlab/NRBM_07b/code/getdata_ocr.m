function [trainData testData] = getdata_ocr(fdata,fold,trainsize)
%function [trainData testData] = getdata_ocr(fdata,fold,trainsize)
%
% load ocr data which was divided in 10 folds

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

	load(fdata);
	trainData	= [];
	testData	= [];
	nfold		= numel(alldata);
	istrain	= zeros(1,nfold*2);
	istrain(fold:(fold+trainsize-1)) = 1;
	istrain = istrain(1:nfold) + istrain((nfold+1):end);
	for f=1:nfold
		if istrain(f)
			if isempty(trainData)
				trainData = alldata{f};
				trainData.X = trainData.X';
			else
				Seq = alldata{f}.Seq + size(trainData.X,2);
				trainData.X	= [trainData.X , alldata{f}.X'];
				trainData.Y	= [trainData.Y ; alldata{f}.Y];
				trainData.Seq	= [trainData.Seq ; Seq];
			end
		else
			if isempty(testData)
				testData = alldata{f};
				testData.X = testData.X';
			else
				Seq = alldata{f}.Seq + size(testData.X,2);
				testData.X	= [testData.X , alldata{f}.X'];
				testData.Y	= [testData.Y ; alldata{f}.Y];
				testData.Seq	= [testData.Seq ; Seq];
			end
		end
	end

