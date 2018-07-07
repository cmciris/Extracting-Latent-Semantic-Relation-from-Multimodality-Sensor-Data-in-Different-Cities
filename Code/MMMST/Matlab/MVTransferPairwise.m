% ************************Multi-view Transfer(our method)*************************
% Created by Ying Wei
% 2015/01/20
% Input: trainS: training exaples from the same distribution
%        trainA: training exaples from the source/different distribution
%        labelS: labels for the TrainS
%        labelA: labels for the TrainA
%        test: test examples
%        groundTruth: labels for the Test
%        N: the number of iterations
%        tag: 0 - our method
%             1 - tradaboost
% Output: Accuracy: the overall accuracy for all the test data
%         PredLabel: the predicted label for all the test data
function [accuracy, predLabel] = MVTransferPairwise(trainS,trainA,labelS,labelA,test,groundTruth,N,tag)
numClass = length(unique([labelA;labelS]));
numTest = size(test,1);
% train pairwise models
pairwise = nchoosek([1:numClass],2);
prob = zeros(numTest,size(pairwise,1));

for k=1:size(pairwise,1)
	idxS = any( bsxfun(@eq, labelS, pairwise(k,:)) , 2 );
	idxA = any( bsxfun(@eq, labelA, pairwise(k,:)) , 2 );	
	idxG = any( bsxfun(@eq, groundTruth, pairwise(k,:)) , 2 );
	labelSK = labelS(idxS,:); trainSK = trainS(idxS,:);
	labelAK = labelA(idxA,:); trainAK = trainA(idxA,:);
% 	groundTruthK = groundTruth(idxG,:); testK = test(idxG,:);
    if(tag == 0)
        [~,prob(:,k)] = MVTrAdaBoost(trainS,trainA,labelSK,labelAK,test,groundTruthK,N,tag);
    else
        [~,prob(:,k)] = TrAdaBoost(trainSK,trainAK,labelSK,labelAK,test,groundTruth,N,tag);
    end
end
% predict the class with the highest probability
predLabel = mode(prob,2);
accuracy = sum(predLabel == groundTruth) ./ numel(groundTruth);    %# accuracy
end