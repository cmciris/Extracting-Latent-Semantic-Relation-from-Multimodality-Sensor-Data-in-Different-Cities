% **************************Experimental Script****************************
% Created by Ying Wei
% 2015/1/13
% 
function [accOur, accB1, accB2, accB3, accB4, accB5] = ExpSCTransferImgDoc(sImg,sDoc,sLabel,tImg,tDoc,tLabel, dictDirectory, ratio, N, rep, sparsity)
    % load the source data from all modalities
%     m1Filename = strcat(sourceDir,'image.txt');
%     m2Filename = strcat(sourceDir,'tag.txt');
%     mlFilename = strcat(sourceDir,'label.txt');
% 
%     sImg = load(m1Filename); 
%     sDoc = load(m2Filename);
%     sLabel = load(mlFilename);
    nSource = size(sLabel,1);
    % load the target data from all modalities
%     n1Filename = strcat(targetDir,'image.txt');
%     n2Filename = strcat(targetDir,'tag.txt');
%     nlFilename = strcat(targetDir,'label.txt');
% 
%     tImg = load(n1Filename);
%     tDoc = load(n2Filename);
%     tLabel = load(nlFilename);
    nTarget = size(tLabel,1);
    sImg = FeaturePreprocess(sImg',0); sImg = sImg';
    sDoc = FeaturePreprocess(sDoc',0); sDoc = sDoc';
    tImg = FeaturePreprocess(tImg',0); tImg = tImg';
    tDoc = FeaturePreprocess(tDoc',0); tDoc = tDoc';
    
%     xRoad = [sRoad; tRoad]; xRoad = FeaturePreprocess(xRoad,0); xRoad(isnan(xRoad))=0; sRoad = xRoad(1:nSource,:); tRoad = xRoad(nSource+1:end,:);
%     xPoi = [sPoi; tPoi]; xPoi = FeaturePreprocess(xPoi,0); xPoi(isnan(xPoi))=0; sPoi = xPoi(1:nSource,:); tPoi = xPoi(nSource+1:end,:);
%     xMet = [sMet; tMet]; xMet = FeaturePreprocess(xMet,0); xMet(isnan(xMet))=0; sMet = xMet(1:nSource,:); tMet = xMet(nSource+1:end,:);
%     sMob = FeaturePreprocess(sMob,0); sMob(isnan(sMob))=0;
    
    % load the dictionary
    dict1Filename = strcat(dictDirectory,'\00.txt');
    dict2Filename = strcat(dictDirectory,'\01.txt');

    D1 = load(dict1Filename); 
    D2 = load(dict2Filename); 

    nCluster = size(D1,1);
    nCluster2 = size(D2,1);

%     nMax = max([nCluster nCluster2]);

    if(nCluster ~= nCluster2)
        nMax = max([D1(:,end-2); D2(:,end-2)])+1;
        fprintf('the dictionary size is not equal!');
        overlap = intersect(D1(:,end-2),D2(:,end-2));
        I = logical(zeros(nMax,1)); I(overlap+1,:)=1;
        I1 = logical(zeros(size(D1,1),1)); idx = arrayfun(@(x) find(D1(:,end-2) == x,1,'first'), overlap ); I1(idx,:)=1;
        I2 = logical(zeros(size(D2,1),1)); idx = arrayfun(@(x) find(D2(:,end-2) == x,1,'first'), overlap ); I2(idx,:)=1;
        D = zeros(nMax,size(D1,2)-3); 
        D(I,:) = D1(I1,1:end-3);  
        if(~isempty(D1(~I1,1:end-3))) 
            tmp = D1(:,end-2)+1;
            D(tmp(~I1,:),:) = D1(~I1,1:end-3); 
        end
        D1 = D;
        D = zeros(nMax,size(D2,2)-3); 
        D(I,:) = D2(I2,1:end-3);  
        if(~isempty(D2(~I2,1:end-3))) 
            tmp = D2(:,end-2)+1;
            D(tmp(~I2,:),:) = D2(~I2,1:end-3);  
        end
        D2 = D;
        
        idx = ~any([D1 D2],2);
        D1(idx,:)=[]; D2(idx,:)=[];
    else
        D1 = D1(:,1:end-3); 
        D2 = D2(:,1:end-3);
    end
    
    nMax=max([size(D1,1) size(D2,1)]);
    
    %hemap
    [h1,h2]=spectralEmbedding(tImg,tDoc,nMax,1);
    tH = [h1 h2]; tH = tH';
    
    % start to sparse coding of the source data
    lambda = [50 5];
%     A = MMSparseCoding(lambda,0,D1,sImg,D2,sDoc);
    sMA = MMSparseCoding(lambda,1,D1,sImg,D2,sDoc,'imgdoc');
    sMA = FeaturePreprocess(sMA',2); sMA = sMA';

%     A = FeaturePreprocess(A',2); A = A';
% %     A(A<0)=0;
%     A1 = A(1:nMax,:);
%     A2 = A(nMax+1:nMax*2,:);

    A1 = MMSparseCoding(lambda(1),0,D1,sImg);
    A2 = MMSparseCoding(lambda(2),0,D2,sDoc);
    A1 = FeaturePreprocess(A1',2); A1 = A1';
    A2 = FeaturePreprocess(A2',2); A2 = A2';

%     A1tmp = A1;
%     A2tmp = A2;
%     A1tmp(A1<A2) = 0;
%     A2tmp(A2<A1) = 0;
%     A1 = A1tmp;
%     A2 = A2tmp;

    % start to sparse coding of the target data
    tA = MMSparseCoding(lambda,1,D1,tImg,D2,tDoc,'imgdoc');
%     t1A = MMSparseCoding(lambda,0,D1,tImg,D2,tDoc);
    tA = FeaturePreprocess(tA',2); tA = tA';
%     t1A = FeaturePreprocess(t1A',2); t1A = t1A';
%     t2A = MMSparseCoding(lambda,0,D3,tMet);
%     tA = [t1A' t2A'];
    
    accB1 = [];
    accB2 = [];
    accB3 = [];
    accB4 = [];
    accB5 = [];
    accOur = [];
    for i = 1:rep
        disp(i);
    % start to randomly split the target data
    nTrain = fix(nTarget * ratio);
    nTest = nTarget - nTrain;
    index = randperm(nTarget,nTrain);
    trIndex = logical(zeros(nTarget,1));
    trIndex(index,:) = 1;
    teIndex = ~trIndex;
    trA = tA(:,trIndex);
    teA = tA(:,teIndex);
    trH = tH(:,trIndex);
    teH = tH(:,teIndex);
    
%     tr1A = t1A(:,trIndex);
%     te1A = t1A(:,teIndex);
    trLabel = tLabel(trIndex,:);
    teLabel = tLabel(teIndex,:);
    % start to multi-view transfer our method
    sA = [A1; A2];
    param.N = N;

%     [accuracy,predLabel] = MVTransfer(trA',sA',trLabel,sLabel,teA',teLabel,param.N,0);
%     accOur = [accOur;accuracy(1)];
    % baseline 1 : multi-feature combined with out transfer
    b1Data = [tImg tDoc];
    b1Train = b1Data(trIndex,:);
    b1Test = b1Data(teIndex,:);
    [acc, pl] = MFeatNoTransferB1(b1Train,trLabel,b1Test,teLabel,0);
    accB1 = [accB1;acc(1)];
    % baseline 2 : sparse codes with out transfer
    [acc, pl] = MFeatNoTransferB1(trA',trLabel,teA',teLabel,0);
    accB2 = [accB2;acc(1)];
    % baseline 3: multi-feature combined Tradaboost
    b3Data = [sImg sDoc];
    [acc,pl] = TrAdaBoost(b1Train,b3Data,trLabel,sLabel,b1Test,teLabel,N,1,1,0);
    accB3 = [accB3;acc(1)];
    % baseline 4: sparse codes Tradaboost
    [acc,pl] = TrAdaBoost(trA',sMA',trLabel,sLabel,teA',teLabel,N,1,1,0);
    accB4 = [accB4; acc(1)];
    % baseline 5: HeMap
    [acc, pl] = MFeatNoTransferB1(trH',trLabel,teH',teLabel,0);
    accB5 = [accB5; acc(1)];
    % our method
    [acc,pl] = MVTrAdaBoost(trA',sA',trLabel,sLabel,teA',teLabel,N,1,1,0);
    accOur = [accOur;acc(1)];
%     disp([num2str(accB2(i)),'-',num2str(accB1(i)),'=',num2str(accB2(i)-accB1(i))]);
%     if(i == 5
    disp([accB1,accB2,accB3,accB4,accB5,accOur]);
%     end
    end
    accOur = mean(accOur);
    accB1 = mean(accB1);
    accB2 = mean(accB2);
    accB3 = mean(accB3);
    accB4 = mean(accB4);
    accB5 = mean(accB5);
end
