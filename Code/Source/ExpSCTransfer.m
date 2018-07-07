% **************************Experimental Script****************************
% Created by Ying Wei
% 2015/1/13
% 
function [accOur] = ExpSCTransfer(sourceDir, targetDir, dictDirectory, timeslot, timeslotS, ratio, N, rep)
	% load the source data from all modalities
	sourceDataDir = strcat(sourceDir, num2str(timeslotS,'%.2d'),'\');
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
	targetDataDir = strcat(targetDir, num2str(timeslot,'%.2d'),'\');
	n1Filename = strcat(targetDataDir,'road.txt');
	n2Filename = strcat(targetDataDir,'poi.txt');
	n3Filename = strcat(targetDataDir,'meterology.txt');
	nlFilename = strcat(targetDataDir,'label.txt');

	tRoad = load(n1Filename);
	tPoi = load(n2Filename);
	tMet = load(n3Filename);
	tLabel = load(nlFilename);
	nTarget = size(tLabel,1);
        
	% load the dictionary
	dict1Filename = strcat(dictDirectory,num2str(timeslotS,'%.2d'),'\00.txt');
	dict2Filename = strcat(dictDirectory,num2str(timeslotS,'%.2d'),'\01.txt');
	dict3Filename = strcat(dictDirectory,num2str(timeslotS,'%.2d'),'\02.txt');
	dict4Filename = strcat(dictDirectory,num2str(timeslotS,'%.2d'),'\03.txt');

	D1 = load(dict1Filename); 
	D2 = load(dict2Filename); 
	D3 = load(dict3Filename); 
	D4 = load(dict4Filename); 
	nCluster = size(D1,1);
    nCluster2 = size(D2,1);
    nCluster3 = size(D3,1);
    nCluster4 = size(D4,1);

    if(nCluster ~= nCluster2 || nCluster ~= nCluster3 || nCluster ~= nCluster4 || nCluster2 ~= nCluster3 || nCluster2 ~= nCluster4 || nCluster3 ~= nCluster4)
         nMax = max([D1(:,end-2); D2(:,end-2); D3(:,end-2); D4(:,end-2)])+1;
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
    
    nMax= max([size(D1,1); size(D2,1);size(D3,1);size(D4,1)]);    

	% start to sparse coding of the source data
	lambda = [0.1 0.1 0.1 0.1];
	A = MMSparseCoding(lambda,0,D1,sRoad,D2,sPoi,D3,sMet,D4,sMob);
    sMA = MMSparseCoding(lambda,1,D1,sRoad,D2,sPoi,D3,sMet,D4,sMob);
    
	A1 = A(1:nMax,:);
	A2 = A(nMax+1:nMax*2,:);
	A3 = A(2*nMax+1:3*nMax,:);
	A4 = A(3*nMax+1:end,:);
    
	% start to sparse coding of the target data
    lambda = [10 10 0.1];
    tA = MMSparseCoding(lambda,1,D1,tRoad,D2,tPoi,D3,tMet);
    
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

        trLabel = tLabel(trIndex,:);
        teLabel = tLabel(teIndex,:);
        % start to multi-view transfer our method
        sA = [A1; A2; A3; A4];
        [acc,pl] = MVTrAdaBoost(trA',sA',trLabel,sLabel,teA',teLabel,N,1,1,1);
        accOur = [accOur;acc(1)];
        disp([timeslot*ones(i,1), accOur]);
    end
    accOur = mean(accOur);
end
