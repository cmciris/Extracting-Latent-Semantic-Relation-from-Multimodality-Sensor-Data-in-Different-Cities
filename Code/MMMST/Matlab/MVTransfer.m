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
function [accuracy, predLabel] = MVTransfer(trainS,trainA,labelS,labelA,test,groundTruth,N,tag)
numClass = length(unique([labelA;labelS]));
numTest = size(test,1);
% train one-against-all models
prob = zeros(numTest,numClass);
for k=1:numClass
	labelSK = ((labelS == k) + 0); labelSK(labelSK == 0) = -1;
	labelAK = ((labelA == k) + 0); labelAK(labelAK == 0) = -1;
	groundTruthK = ((groundTruth == k) + 0); groundTruthK(groundTruthK == 0) = -1;
    if(tag == 0)
        [~,prob(:,k)] = MVTrAdaBoost(trainS,trainA,labelSK,labelAK,test,groundTruthK,N,0);
    else
        [prob(:,k),~] = TrAdaBoost(trainS,trainA,labelSK,labelAK,test,groundTruthK,N,0);
%         ClassTreeEns = fitensemble(trainS,labelSK,'AdaBoostM1',100,'Tree');
%         prob(:,k) = predict(ClassTreeEns,test);
    end
end
% predict the class with the highest probability
[~,predLabel] = max(prob,[],2);
accuracy = sum(predLabel == groundTruth) / numel(groundTruth);    %# accuracy
end