% **************************Experimental Script****************************
% Created by Ying Wei
% 2015/1/13
% 
function [accB1, accB2] = TestSparsecoding(sourceDir, targetDir, dictDirectory, timeslot, timeslotD, ratio, N, rep, sparsity)
	% load the source data from all modalities
	sourceDataDir = strcat(sourceDir,'labeledSubset4\',num2str(timeslot,'%.2d'),'\');
	m1Filename = strcat(sourceDataDir,'road.txt');
	m2Filename = strcat(sourceDataDir,'poi.txt');
	m3Filename = strcat(sourceDataDir,'meterology.txt');
	m4Filename = strcat(sourceDataDir, 'mobility.txt');
	mlFilename = strcat(sourceDataDir,'label.txt');

	sRoad = load(m1Filename); 
	sPoi = load(m2Filename);
	sMet = load(m3Filename);
	sMob = load(m4Filename);
	sLabel = load(mlFilename);
	nSource = size(sLabel,1);
	% load the target data from all modalities
	targetDataDir = strcat(targetDir,'labeledSubset4\',num2str(timeslot,'%.2d'),'\');
	n1Filename = strcat(targetDataDir,'road.txt');
	n2Filename = strcat(targetDataDir,'poi.txt');
	n3Filename = strcat(targetDataDir,'meterology.txt');
	nlFilename = strcat(targetDataDir,'label.txt');

	tRoad = load(n1Filename);
	tPoi = load(n2Filename);
	tMet = load(n3Filename);
	tLabel = load(nlFilename);
	nTarget = size(tLabel,1);
    
%     xRoad = [sRoad; tRoad]; xRoad = FeaturePreprocess(xRoad,0); xRoad(isnan(xRoad))=0; sRoad = xRoad(1:nSource,:); tRoad = xRoad(nSource+1:end,:);
%     xPoi = [sPoi; tPoi]; xPoi = FeaturePreprocess(xPoi,0); xPoi(isnan(xPoi))=0; sPoi = xPoi(1:nSource,:); tPoi = xPoi(nSource+1:end,:);
%     xMet = [sMet; tMet]; xMet = FeaturePreprocess(xMet,0); xMet(isnan(xMet))=0; sMet = xMet(1:nSource,:); tMet = xMet(nSource+1:end,:);
%     sMob = FeaturePreprocess(sMob,0); sMob(isnan(sMob))=0;

% D11 = []; D22 = []; D33 = []; D44 = [];
% for timeslotD = 0:23
	% load the dictionary
	dict1Filename = strcat(dictDirectory,num2str(timeslotD,'%.2d'),'\00.txt');
	dict2Filename = strcat(dictDirectory,num2str(timeslotD,'%.2d'),'\01.txt');
	dict3Filename = strcat(dictDirectory,num2str(timeslotD,'%.2d'),'\02.txt');
	dict4Filename = strcat(dictDirectory,num2str(timeslotD,'%.2d'),'\03.txt');

	D1 = load(dict1Filename); 
	D2 = load(dict2Filename); 
	D3 = load(dict3Filename); 
	D4 = load(dict4Filename); 
	nCluster = size(D1,1);
    nCluster2 = size(D2,1);
    nCluster3 = size(D3,1);
    nCluster4 = size(D4,1);
    nMax = max([D1(:,end-2); D2(:,end-2); D3(:,end-2); D4(:,end-2)])+1;

    if(nCluster ~= nCluster2 || nCluster ~= nCluster3 || nCluster ~= nCluster4 || nCluster2 ~= nCluster3 || nCluster2 ~= nCluster4 || nCluster3 ~= nCluster4)
        fprintf('the dictionary size is not equal!');
        overlap = intersect(D1(:,end-2),intersect(D2(:,end-2),intersect(D3(:,end-2),D4(:,end-2))));
        I = logical(zeros(nMax,1)); I(overlap+1,:)=1;
        I1 = logical(zeros(size(D1,1),1)); idx = arrayfun(@(x) find(D1(:,end-2) == x,1,'first'), overlap ); I1(idx,:)=1;
        I2 = logical(zeros(size(D2,1),1)); idx = arrayfun(@(x) find(D2(:,end-2) == x,1,'first'), overlap ); I2(idx,:)=1;
        I3 = logical(zeros(size(D3,1),1)); idx = arrayfun(@(x) find(D3(:,end-2) == x,1,'first'), overlap ); I3(idx,:)=1;
        I4 = logical(zeros(size(D4,1),1)); idx = arrayfun(@(x) find(D4(:,end-2) == x,1,'first'), overlap ); I4(idx,:)=1;
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
		D = zeros(nMax,size(D3,2)-3); 
        D(I,:) = D3(I3,1:end-3);  
        if(~isempty(D3(~I3,1:end-3))) 
            tmp = D3(:,end-2)+1;
            D(tmp(~I3,:),:) = D3(~I3,1:end-3); 
        end
        D3 = D;
		D = zeros(nMax,size(D4,2)-3); 
        D(I,:) = D4(I4,1:end-3);  
        if(~isempty(D4(~I4,1:end-3))) 
            tmp = D4(:,end-2)+1;
            D(tmp(~I4,:),:) = D4(~I4,1:end-3);  
        end
        D4 = D;
        idx = ~any([D1 D2 D3 D4],2);
        D1(idx,:)=[]; D2(idx,:)=[]; D3(idx,:)=[]; D4(idx,:)=[];
    else
        D1 = D1(:,1:end-3); 
        D2 = D2(:,1:end-3);
        D3 = D3(:,1:end-3);
        D4 = D4(:,1:end-3);        
    end
%     D11 = [D11; D1];
%     D22 = [D22; D2];
%     D33 = [D33; D3];
%     D44 = [D44; D4];
% 
% end
%     if(nCluster ~= nCluster2)
%         overlap = intersect(D1(:,end-2),D2(:,end-2));
%         I1 = arrayfun(@(x) find(D1(:,end-2) == x,1,'first'), overlap );
%         I2 = arrayfun(@(x) find(D2(:,end-2) == x,1,'first'), overlap );
%         D1 = D1(I1,1:end-3); 
%         D2 = D2(I2,1:end-3);
%     else
%         D1 = D1(:,1:end-3); 
%         D2 = D2(:,1:end-3);
%     end
%     if(nCluster3 ~= nCluster4)
%         overlap = intersect(D3(:,end-2),D4(:,end-2));
%         I3 = arrayfun(@(x) find(D3(:,end-2) == x,1,'first'), overlap );
%         I4 = arrayfun(@(x) find(D4(:,end-2) == x,1,'first'), overlap ); 
%         D3 = D3(I3,1:end-3);
%         D4 = D4(I4,1:end-3);
%     else
%         D3 = D3(:,1:end-3);
%         D4 = D4(:,1:end-3);         
%     end
    nMax=max([size(D1,1) size(D2,1) size(D3,1) size(D4, 1)]);

	% start to sparse coding of the source data
	lambda = [1 1 sparsity sparsity];
	A = MMSparseCoding(lambda,0,D1,sRoad,D2,sPoi,D3,sMet,D4,sMob);
    sMA = MMSparseCoding(lambda,1,D1,sRoad,D2,sPoi,D3,sMet,D4,sMob);
    
%     A = FeaturePreprocess(A',2); A = A';
	A1 = A(1:nMax,:);
	A2 = A(nMax+1:nMax*2,:);
	A3 = A(2*nMax+1:3*nMax,:);
	A4 = A(3*nMax+1:end,:);
    
	% start to sparse coding of the target data
    lambda = [10 10 0.1];
    tA = MMSparseCoding(lambda,1,D1,tRoad,D2,tPoi,D3,tMet);
%     t1A = MMSparseCoding(lambda,1,D1,tRoad,D2,tPoi);
%     t2A = MMSparseCoding(lambda,0,D3(:,4:5),tMet(:,4:5));
%     sMA = FeaturePreprocess(sMA',1); sMA = sMA';
%     tA = FeaturePreprocess(tA',0); tA = tA';
%     tA = [t1A' t2A'];
%     t2A = tA(any(tA,2),:);
    
    accB1 = [];
    accB2 = [];
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
%     tr2A = t2A(:,trIndex);
% 	te2A = t2A(:,teIndex);
% 	trA = tA(trIndex,:);
% 	teA = tA(teIndex,:);
	trLabel = tLabel(trIndex,:);
	teLabel = tLabel(teIndex,:);
    % start to multi-view transfer our method
    sA = [A1; A2; A3; A4];
    param.N = N;
%     [accuracy,predLabel] = MVTransfer(trA',sA',trLabel,sLabel,teA',teLabel,param.N,0);
%     accOur = [accOur;accuracy(1)];
    % baseline 1 : multi-feature combined with out transfer
    b1Data = [tRoad tPoi tMet];
    b1Train = b1Data(trIndex,:);
    b1Test = b1Data(teIndex,:);
    [acc, pl] = MFeatNoTransferB1(b1Train,trLabel,b1Test,teLabel,1);
    accB1 = [accB1;acc(1)];
    % baseline 2 : sparse codes with out transfer
    [acc, pl] = MFeatNoTransferB1(trA',trLabel,teA',teLabel,1);
    accB2 = [accB2;acc(1)];
    end
    accB1 = mean(accB1);
    accB2 = mean(accB2);
end
