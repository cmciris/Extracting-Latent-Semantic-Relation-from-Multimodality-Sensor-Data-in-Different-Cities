function [accuracy,predLabels] = CoTraining(train, label, test, groundTruth, noDim1, topK, varargin)
    predLabels = [];
   	index = [];
   	allIndex = [1:size(test,1)]';
    count = 0;
    if(nargin >= 7)
        trainWeight = varargin{1,1};
        nT = varargin{1,2};
        weight = [trainWeight; ones(size(test(nT:end,:),1),1)];
    end
	while(~isempty(test))
%         count = count + 1;
%         disp(count);
		% split the training data into two sources
		train1 = train(:,1:noDim1);
		train2 = train(:,noDim1+1:end);
		test1 = test(:,1:noDim1);
		test2 = test(:,noDim1+1:end);


		if(nargin <7)
			% train a classifier using the first source
			tree1 = fitctree(train1,label);
			% train a classifier using the second source
			tree2 = fitctree(train2,label);
		else
			% train a classifier using the first source
			tree1 = fitctree(train1,label,'Weights', trainWeight);
			% train a classifier using the second source
			tree2 = fitctree(train2,label,'Weights', trainWeight);
		end
			

		% predict the result 
		[label1, score1] = predict(tree1, test1);
		[label2, score2] = predict(tree2, test2);

		% pick the top k confidently predicted test instances and add them to the training pool
		[~, I1] = sort(max(score1,[],2),'descend');
		[~, I2]	= sort(max(score2,[],2),'descend');
        if(size(I1,1)<topK)
            topK = size(I1,1);
        end
		remove = union(I1(1:topK),I2(1:topK));
		predLabel = zeros(size(test,1),1);
        tmpLabel = zeros(size(test,1),2);
		tmpLabel(I1(1:topK),1) = label1(I1(1:topK),:);
		tmpLabel(I2(1:topK),2) = label2(I2(1:topK),:);
		[~,I] = sort([max(score1,[],2) max(score2,[],2)],2);
        for i =1:size(test,1)
            if(tmpLabel(i,I(i,2))>0)
                predLabel(i) = tmpLabel(i,I(i,2));
            else
                predLabel(i) = tmpLabel(i,I(i,1));
            end
        end
		predLabel(predLabel==0,:)=[]; 
		train = [train; test(remove,:)];
		if(nargin >= 7)
			trainWeight = [trainWeight; weight(remove,:)];
			weight(remove,:)=[];
		end
		label = [label; predLabel];
        predLabels = [predLabels;predLabel];
        test(remove,:)=[];
		index = [index; allIndex(remove,:)];
		allIndex(remove,:)=[];
    end
    [~,sorted]=sort(index);
    predLabels = predLabels(sorted,:);
    tmp = predLabels == groundTruth;
    accuracy = sum(tmp)/numel(tmp);
end