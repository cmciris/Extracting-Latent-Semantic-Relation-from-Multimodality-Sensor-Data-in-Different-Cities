% ****************Baseline 1: Multi-Feature Combination Without Transfer************
% Created by Ying Wei
% 2015/1/19
% Input: train: training data inã€?he target domain with each instance naively combining all the features
% 	     trLabelï¼?training label in the target domain
%	 	 test: test data in the target domain
%		 teLabel: test label in the target domain
% Output: accuracy: the classification of the accuracy on the test data
% 		  predLabel: the predicted label for the test instances 
function  [accuracy, predLabel] = MFeatNoTransferB1(train, trLabel, test, teLabel, tag)
	 addpath libsvm-3.20\matlab
     addpath liblinear-1.96\liblinear-1.96\matlab
	if(tag == 0)
%         tic
% 		model = svmtrain(trLabel, train, '-t 0 -c 10 -q');
% 		[predLabel, accuracy, dec_values] = svmpredict(teLabel, test, model);
%         toc
		tree = fitctree(train,trLabel);
		% view(tree,'Mode','graph')
		predLabel = predict(tree,test);
%         predLabel = round(predLabel);
		error = (predLabel ~= teLabel) + 0;
		accuracy = 1 - sum(error) / numel(error);
	else
		tree = fitctree(train,trLabel);
		% view(tree,'Mode','graph')
		predLabel = predict(tree,test);
		error = (predLabel ~= teLabel) + 0;
		accuracy = 1 - sum(error) / numel(error);
	end
end