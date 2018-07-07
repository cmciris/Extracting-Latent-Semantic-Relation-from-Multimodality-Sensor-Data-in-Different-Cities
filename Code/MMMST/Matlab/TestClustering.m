function str = TestClustering(targetDir,dictDirectory,timeslot,lambda)
	% load the target data from all modalities
	targetDataDir = strcat(targetDir,'labeledSubset\',num2str(timeslot,'%.2d'),'\');
	n1Filename = strcat(targetDataDir,'road.txt');
	n2Filename = strcat(targetDataDir,'poi.txt');
	n3Filename = strcat(targetDataDir,'meterology.txt');
%     n4Filename = strcat(targetDataDir,'mobility.txt');
	nlFilename = strcat(targetDataDir,'label.txt');

	tRoad = load(n1Filename);
	tPoi = load(n2Filename);
	tMet = load(n3Filename);
%     tMob = load(n4Filename);
	tLabel = load(nlFilename);
	nTarget = size(tLabel,1);
	b1Data = [tRoad tPoi tMet];

	% load the dictionary
	dict1Filename = strcat(dictDirectory,num2str(timeslot,'%.2d'),'\00.txt');
	dict2Filename = strcat(dictDirectory,num2str(timeslot,'%.2d'),'\01.txt');
	dict3Filename = strcat(dictDirectory,num2str(timeslot,'%.2d'),'\02.txt');
	dict4Filename = strcat(dictDirectory,num2str(timeslot,'%.2d'),'\03.txt');

	D1 = load(dict1Filename); 
	D2 = load(dict2Filename); 
	D3 = load(dict3Filename); 
	D4 = load(dict4Filename); 

% 	nCluster = size(D1,1);
%     nCluster2 = size(D2,1);
%     nCluster3 = size(D3,1);
%     nCluster4 = size(D4,1);
%     nMax = max([D1(:,end-2); D2(:,end-2); D3(:,end-2); D4(:,end-2)])+1;
% 
%     if(nCluster ~= nCluster2 || nCluster ~= nCluster3 || nCluster ~= nCluster4 || nCluster2 ~= nCluster3 || nCluster2 ~= nCluster4 || nCluster3 ~= nCluster4)
%         fprintf('the dictionary size is not equal!');
%         overlap = intersect(D1(:,end-2),intersect(D2(:,end-2),intersect(D3(:,end-2),D4(:,end-2))));
%         I = logical(zeros(nMax,1)); I(overlap+1,:)=1;
%         I1 = logical(zeros(size(D1,1),1)); idx = arrayfun(@(x) find(D1(:,end-2) == x,1,'first'), overlap ); I1(idx,:)=1;
%         I2 = logical(zeros(size(D2,1),1)); idx = arrayfun(@(x) find(D2(:,end-2) == x,1,'first'), overlap ); I2(idx,:)=1;
%         I3 = logical(zeros(size(D3,1),1)); idx = arrayfun(@(x) find(D3(:,end-2) == x,1,'first'), overlap ); I3(idx,:)=1;
%         I4 = logical(zeros(size(D4,1),1)); idx = arrayfun(@(x) find(D4(:,end-2) == x,1,'first'), overlap ); I4(idx,:)=1;
%         D = zeros(nMax,size(D1,2)-3); 
%         D(I,:) = D1(I1,1:end-3);  
%         if(~isempty(D1(~I1,1:end-3))) 
%             tmp = D1(:,end-2)+1;
%             D(tmp(~I1,:),:) = D1(~I1,1:end-3); 
%         end
%         D1 = D;
% 		D = zeros(nMax,size(D2,2)-3); 
%         D(I,:) = D2(I2,1:end-3);  
%         if(~isempty(D2(~I2,1:end-3))) 
%             tmp = D2(:,end-2)+1;
%             D(tmp(~I2,:),:) = D2(~I2,1:end-3);  
%         end
%         D2 = D;
% 		D = zeros(nMax,size(D3,2)-3); 
%         D(I,:) = D3(I3,1:end-3);  
%         if(~isempty(D3(~I3,1:end-3))) 
%             tmp = D3(:,end-2)+1;
%             D(tmp(~I3,:),:) = D3(~I3,1:end-3); 
%         end
%         D3 = D;
% 		D = zeros(nMax,size(D4,2)-3); 
%         D(I,:) = D4(I4,1:end-3);  
%         if(~isempty(D4(~I4,1:end-3))) 
%             tmp = D4(:,end-2)+1;
%             D(tmp(~I4,:),:) = D4(~I4,1:end-3);  
%         end
%         D4 = D;
%         idx = ~any([D1 D2 D3 D4],2);
%         D1(idx,:)=[]; D2(idx,:)=[]; D3(idx,:)=[]; D4(idx,:)=[];
%     else
%         D1 = D1(:,1:end-3); 
%         D2 = D2(:,1:end-3);
%         D3 = D3(:,1:end-3);
%         D4 = D4(:,1:end-3);        
%     end

	% start to sparse coding of the target data
% 	tA = MMSparseCoding(lambda,1,D1,tRoad,D2,tPoi,D3,tMet);
	tA = MMSparseCoding(lambda,0,D3(:,1:end-3),tMet);

    tA(tA == -1e+12)=0;
	rep = 30;

	accAll1 = []; accAll2 = [];
	for ratio = 0.1:0.1:1
		% start to randomly split the target data
		for r = 1:rep
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

		accuracy1 = [];
		accuracy2 = [];
		[acc, pl] = MFeatNoTransferB1(b1Data(trIndex,:),trLabel,b1Data,tLabel);
		accuracy1 = [accuracy1 acc];
		[acc, pl] = MFeatNoTransferB1(trA',trLabel,tA',tLabel);
		accuracy2 = [accuracy2 acc];
		end
		accAll1 = [accAll1 mean(accuracy1)];
		accAll2 = [accAll2 mean(accuracy2)];
    end
	plot(accAll1,'b');
    hold on;
	plot(accAll2);
    hold on;
end