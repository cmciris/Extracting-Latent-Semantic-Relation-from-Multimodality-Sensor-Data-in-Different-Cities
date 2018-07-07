% % ******************************Tradaboost********************************
% Created by Ying Wei
% 2015/1/5
% Input: trainS: training exaples from the same distribution
%        trainA: training exaples from the source/different distribution
%        labelS: labels for the TrainS
%        labelA: labels for the TrainA
%        test: test examples
%        groundTruth: labels for the Test
%        N: the number of iterations
% Output: Accuracy: the overall accuracy for all the test data
%         PredLabel: the predicted label for all the test data
function [accAll, predictLabels] = TrAdaBoost(trainS,trainA,labelS,labelA,test,groundTruth,N,initWeightA,initWeightS,tag)
addpath libsvm-weights-3.18\matlab

trainData = [trainA;trainS];
trainLabel = [labelA;labelS];

numClass = length(unique(labelA));

[nTrainS,dimS] = size(trainS);
[nTrainA,dimA] = size(trainA);
[nTest,dimT] = size(test);

accAll = [];
%initial the weights of the instances
for alpha = initWeightA%[0.001 0.01 0.1 1 10 100 1000]
    for beta = initWeightS%[0.001 0.01 0.1 1 10 100 1000]
weight = [alpha*ones(nTrainA,1); beta*ones(nTrainS,1)];
beta = 1.0 / (1+sqrt(2*log(nTrainA/N)));% * 0.1;
betaT = zeros(1,N);
resultLabels = ones(nTest,N);
predictLabels = ones(nTest,N);

errorTest = [];
for t=1:N
%     disp([t]);
    %set the weight
    p = weight/sum(weight);
    %call the weak learner svm
    if(tag == 0)
%         tic
%         model = svmtrain(p, trainLabel, trainData, '-t 0 -c 10 -q');
%         [realLabel, accuracy, predict_label] = svmpredict([trainLabel;groundTruth], [trainData;test], model);
%         toc
        tree = fitctree(trainData,trainLabel,'Weights',p);
       	realLabel = predict(tree,[trainData;test]);
%         realLabel = round(realLabel);
    else
        tree = fitctree(trainData,trainLabel,'Weights',p);
        realLabel = predict(tree,[trainData;test]);
    end
%     [realLabel, accuracy, predict_label] = svmpredict([groundTruth], [test], model);
%     if(tag)
%         predict_label = realLabel;
%     else
%         predict_label = predict_label * (2 * model.Label == 1 - 1);
%     end
    resultLabels(:,t) = realLabel(nTrainA + nTrainS+1:end,:);
    predictLabels(:,t) = realLabel(nTrainA+nTrainS+1:end);
    %update betaT
    error = realLabel((nTrainA+1):(nTrainA+nTrainS))~=trainLabel((nTrainA+1):(nTrainA+nTrainS))+0;
    error2 = realLabel((nTrainA+nTrainS+1):end)==groundTruth(1:end)+0;
    betaT(t) = sum(weight((nTrainA+1):(nTrainA+nTrainS)).* error)/sum(weight((nTrainA+1):(nTrainS+nTrainA)));
    errorTest = [errorTest;sum(error2)/numel(error2)]; 
    if(betaT(t)>0.5)
        betaT(t)=0.5;
    end
    if(betaT(t)==0)
        break;
%         betaT(t)=0.01;
    end
    betaT(t) = betaT(t)/(1-betaT(t));
    
    if(betaT(t) >= 1 && t > 1)
        break;
    end

    %upadate the new weight vector
    weight(1:nTrainA)=weight(1:nTrainA).*(beta.^(((realLabel(1:nTrainA)~=trainLabel(1:nTrainA)+0))));
%     weight(1:nTrainA)=weight(1:nTrainA).*(beta.^(abs(realLabel(1:nTrainA)-trainLabel(1:nTrainA))/2));
    weight(nTrainA+1:end)=weight(nTrainA+1:end).*(betaT(t).^(-error));   
%     weight(nTrainA+1:end)=weight(nTrainA+1:end)/betaT(t)/5;
    % tackling the class imbalance problem
%     pIndex = trainLabel == 1;
%     nIndex = ~pIndex;
%     pWeight = sum(weight(pIndex,:));
%     nWeight = sum(weight(nIndex,:));
%     weight(pIndex,:) = weight(pIndex,:) * (pWeight + nWeight) / pWeight / 2; 
%     weight(nIndex,:) = weight(nIndex,:) * (pWeight + nWeight) / nWeight / 2;

%     allWeight = sum(weight);
%     for i =1:numClass
%         iIndex = trainLabel == i;
%         iWeight = sum(weight(iIndex,:));
%         weight(iIndex,:) = weight(iIndex,:) * allWeight /iWeight/numClass;
%     end
end

% disp(errorTest);

if t>1
    finalN = t-1;
else
    finalN = 1;
end
%output the final result of the test data
probEstimates = resultLabels;
% for t = 1:finalN
%     probEstimates = [probEstimates resultLabels(:,t)];
% end
% [~,idx]=sort(betaT(1:finalN));
% idx = idx(1:ceil(finalN/2));
idx =1:finalN;
probEstimates = probEstimates(:,idx); 
% probEstimates = mean(probEstimates,2);
betaT = 0.5*log(1./betaT);

%% solution 1: weighted average of the label and round
% probEstimates = sum(repmat(betaT(idx),nTest,1) .* probEstimates,2) / sum(betaT(idx)); 
% predictLabels = round(probEstimates);
%% solution 2: only consider top-k confident betaT and output the vote label
[~,index] = sort(betaT(idx),'descend');

acc = [];
for k=10:10:N
    if(k>finalN)
       k = finalN;
    end
    tmp = probEstimates(:,index(1:k));
    predictLabels = mode(tmp,2);
    acc = [acc; sum((predictLabels==groundTruth) + 0)/nTest];
end
acc = max(acc);
accAll = [accAll; acc];
    end
end
end