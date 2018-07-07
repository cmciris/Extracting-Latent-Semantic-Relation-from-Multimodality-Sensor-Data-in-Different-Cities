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

% trainA(trainA <0) = 0;

[nTrainS,~]=size(trainS);
[nTest,dim]=size(test);
[nTrainA,dimA]=size(trainA);
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
    %set the weight
    weight = [weightA;repmat(weightS,1,numMod)];
    for i=1:numMod
        weight(:,i)=weight(:,i)/sum(weight(:,i));  %normalization
    end
    tmpLabel = ones(nTrainA+nTrainS,numMod);
    predictLabel = ones(nTrainA+nTrainS,numMod);
    
    for i =1:numMod
        %call the weak learner
        trainLabel = [labelA;labelS];
        trainData = [trainA(:,(i-1)*dim+1:i*dim);trainS];
        if(tag == 0)
            tree = fitctree(trainData,trainLabel,'Weights',weight(:,i));
            [realLabel, score] = predict(tree,[trainData;test]);  
            [tmp,~] =sort(score(nTrainA+nTrainS+1:end,:),2,'descend');
            scores(:,i,t) = tmp(:,1);
        else
            tree = fitctree(trainData,trainLabel,'Weights',weight(:,i));
            [realLabel,score] = predict(tree,[trainData;test]);
            [tmp,~] =sort(score(nTrainA+nTrainS+1:end,:),2,'descend');
            scores(:,i,t) = tmp(:,1);
        end
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
%         break;
        betaT(t)=0.01;
    end
    betaT(t) = betaT(t)/(1-betaT(t));

    if(betaT(t) >= 1 && t > 1)
        break;
    end

    %update the new weight vector
    for i=1:numMod
        weightA(:,i)=weightA(:,i).*(beta.^(abs(tmpLabel(1:nTrainA,i)~=labelA+0)));
    end
    weightS = weightS.*(betaT(t).^(-abs(error)));
end


if t>1
    finalN = t-1;
else
    finalN = 1;
end
%output the final result of the test data
probEstimates=[];
for t = 1:finalN
    if(tag == 0)
        [~,tmpI] = sort(scores(:,:,t),2,'descend');
        tmp = zeros(nTest,2);
        tmp(tmpI(:,1)==1,1) = 1;
        tmp(tmpI(:,1)==2,2) = 1;
        probEstimates = [probEstimates sum(tmp .* resultLabels(:,:,t),2)];
    else
        probEstimates = [probEstimates mode(resultLabels(:,:,t),2)];
    end
end

idx =1:finalN;
probEstimates = probEstimates(:,idx); 
betaT = 0.5*log(1./betaT);

%% solution: only consider top-k confident betaT and output the vote label
[~,index] = sort(betaT(idx),'descend');

acc = [];
for k=10:10:N
    if(k>finalN)
       k = finalN;
    end
    tmp = probEstimates(:,index(1:k));
    predictLabels = mode(tmp,2);
    acc = [acc;sum((predictLabels==groundTruth) + 0)/nTest];
end
acc = max(acc);
end