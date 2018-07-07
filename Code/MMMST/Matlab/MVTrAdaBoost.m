% *****************************MV-Tradaboost*****************************
% Created by Ying Wei
% 2015/1/6
% Input: trainS: training exaples from the same distribution
%        trainA: training exaples from the source/different distribution
%        labelS: labels for the TrainS
%        labelA: labels for the TrainA
%        test: test examples
%        groundTruth: labels for the Test
%        N: the number of iterations
% Output: Accuracy: the overall accuracy for all the test data
%         PredLabel: the predicted label for all the test data
function [acc,predictLabels] = MVTrAdaBoost(trainS,trainA,labelS,labelA,test,groundTruth,N,initWeightA,initWeightS,tag)
addpath libsvm-weights-3.18\matlab

% trainA(trainA <0) = 0;

[nTrainS,dim]=size(trainS);
[nTest,dim]=size(test);
[nTrainA,dimA]=size(trainA);
numClass = size(unique(labelA),1);
numMod = dimA / dim;
%initial the weights
weightA = ones(nTrainA,1)* (initWeightA * ones(1,numMod));
weightS = initWeightS * ones(nTrainS,1);
beta = 1.0 / (1+sqrt(2*log(nTrainA/N)));% * 0.5;
betaT = zeros(1,N);
resultLabels = ones(nTest,numMod,N);

error2Test = [];
scores = zeros(nTest,numMod,N);
%start to iterate and boost
for t=1:N
%     disp([t]);
    % tackling the class imbalance problem
%     psIndex = labelS == 1; paIndex = labelA == 1;
%     nsIndex = ~psIndex; naIndex = ~paIndex;
%     for i =1:numMod
%         pWeight = sum(weightA(paIndex,i)); pWeight = pWeight + sum(weightS(psIndex,:));
%         nWeight = sum(weightA(naIndex,i)); nWeight = nWeight + sum(weightS(nsIndex,:));
%         weightA(paIndex,i) = weightA(paIndex,i) * (pWeight + nWeight) / pWeight / 2; 
%         weightA(naIndex,i) = weightA(naIndex,i) * (pWeight + nWeight) / nWeight / 2;
%         weightS(psIndex,:) = weightS(psIndex,:) * (pWeight + nWeight) / pWeight / 2; 
%         weightS(nsIndex,:) = weightS(nsIndex,:) * (pWeight + nWeight) / nWeight / 2;
%     end
    %set the weight
    weight = [weightA;repmat(weightS,1,numMod)];
    for i=1:numMod
        weight(:,i)=weight(:,i)/sum(weight(:,i));  %normalization
    end
    tmpLabel = ones(nTrainA+nTrainS,numMod);
    predictLabel = ones(nTrainA+nTrainS,numMod);
    
    for i =1:numMod
        %call the weak learner svm
        trainLabel = [labelA;labelS];
        trainData = [trainA(:,(i-1)*dim+1:i*dim);trainS];
        if(tag == 0)
%            	tic;
%             model = svmtrain(weight(:,i), trainLabel, trainData, '-t 0 -c 1 -b 1 -q');
%             [realLabel, accuracy, predictlabel] = svmpredict([trainLabel;groundTruth], [trainData;test], model, '-b 1');
            tree = fitctree(trainData,trainLabel,'Weights',weight(:,i));
            [realLabel, score] = predict(tree,[trainData;test]);  
%             realLabel = round(realLabel);
%             toc;
            [tmp,~] =sort(score(nTrainA+nTrainS+1:end,:),2,'descend');
            scores(:,i,t) = tmp(:,1);
        else
            tree = fitctree(trainData,trainLabel,'Weights',weight(:,i));
            [realLabel,score] = predict(tree,[trainData;test]);
            [tmp,~] =sort(score(nTrainA+nTrainS+1:end,:),2,'descend');
            scores(:,i,t) = tmp(:,1);
        end
%         [~, acc, ~] = svmpredict(groundTruth, test, model,'-b 1');
%         predictLabel(:,i) = realLabel(1:nTrainA+nTrainS,:);
        %this place may need choose the right label's corresponding
        %prob-estimates
        tmpLabel(:,i)= realLabel(1:nTrainA+nTrainS,:);
        resultLabels(:,i,t)=realLabel(nTrainA+nTrainS+1:end,:);
    end
    %calculate the error rate and agreement
    error = max(tmpLabel((nTrainA+1):(nTrainA+nTrainS),:)~=repmat(labelS,1,numMod)+0,[],2);
    error2 = resultLabels(:,:,t)==repmat(groundTruth,1,numMod)+0;
    error2Test = [error2Test; sum(error2) ./ 778];
    betaT(t) = sum(weightS.*abs(error)/sum(weightS));
    agreement = zeros(nTrainS,1);
    for i=1:numMod-1
        for j=i+1:numMod
            agreement = agreement + abs(tmpLabel((nTrainA+1):(nTrainA+nTrainS),i)~=tmpLabel((nTrainA+1):(nTrainA+nTrainS),j)+0);
        end
    end
    agreement = sum(agreement)/((nTrainS) * nchoosek(numMod,2));
    betaT(t) = betaT(t)*(1-agreement);
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

    
    % betaT(t) = betaT(t)/(1-betaT(t)/(numClass-1));
    %update the new weight vector
    for i=1:numMod
        weightA(:,i)=weightA(:,i).*(beta.^(abs(tmpLabel(1:nTrainA,i)~=labelA+0)));
    end
    weightS = weightS.*(betaT(t).^(-abs(error)));
    
%     for i=1:numMod
%         allWeight = sum([weightA(:,i); weightS]);
%         tmp = [weightA(:,i); weightS];
%         for j =1:numClass
%             jIndex = trainLabel == j;
%             jWeight = sum(tmp(jIndex,:));
%             tmp(jIndex,:) = tmp(jIndex,:) * allWeight /jWeight/numClass;
%         end
%         weightA(:,i)= tmp(1:nTrainA,:);
%         weightS = tmp(nTrainA+1:end,:);
%    end
end


if t>1
    finalN = t-1;
else
    finalN = 1;
end
%output the final result of the test data
% probEstimates = mode(resultLabels,2);
% probEstimates = reshape(probEstimates,[],t);
probEstimates=[];
for t = 1:finalN
    if(tag == 0)
        [~,tmpI] = sort(scores(:,:,t),2,'descend');
        tmp = zeros(nTest,2);
        tmp(tmpI(:,1)==1,1) = 1;
        tmp(tmpI(:,1)==2,2) = 1;
        probEstimates = [probEstimates sum(tmp .* resultLabels(:,:,t),2)];
    else
%         [~,tmpI] = sort(scores(:,:,t),2,'descend');
%         tmp = zeros(nTest,4);
%         tmp(tmpI(:,1)==1,1) = 1;
%         tmp(tmpI(:,1)==2,2) = 1;
%         tmp(tmpI(:,1)==3,3) = 1;
%         tmp(tmpI(:,1)==4,4) = 1;
%         probEstimates = [probEstimates sum(tmp .* resultLabels(:,:,t),2)];
        probEstimates = [probEstimates mode(resultLabels(:,:,t),2)];
    end
end
% probEstimates = mean(probEstimates,2);
% acc = sum(probEstimates > 0.5 == groundTruth + 0) / nTest;

% idx =ceil(finalN/2):finalN;
idx =1:finalN;
probEstimates = probEstimates(:,idx); 
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
%     if(k>ceil(finalN/2))
%         k = ceil(finalN/2);
%     end
    tmp = probEstimates(:,index(1:k));
    predictLabels = mode(tmp,2);
    acc = [acc;sum((predictLabels==groundTruth) + 0)/nTest];
end
acc = max(acc);
% temp = ones(nTest,1);
% for t=ceil(N/2):N
%     temp = temp .* (betaT(t).^(-max(resultLabels(:,:,t),[],2)));
%     for i = 1:numClass
%         classes(i)=classes(i) * betaT(t)^(-i);
%     end
% end

% [B,index] = min(abs(repmat(temp,1,numClass)-repmat(classes',nTest,1)),[],2);

% accuracy = sum((index - groundTruth) == 0 + 0)/nTest;
% predLabel = index;
end