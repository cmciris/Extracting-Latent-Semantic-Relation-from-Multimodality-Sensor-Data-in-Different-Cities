function result = AdjustParameters()
	result = [];
    i=1;
    sourceDir = 'D:\Dropbox\Dropbox\Project\MSRA\ImgDocResult\Features\Source\car_vs_dog\';
    targetDir = 'D:\Dropbox\Dropbox\Project\MSRA\ImgDocResult\Features\Target\car_vs_dog\';
    m1Filename = strcat(sourceDir,'image.txt');
    m2Filename = strcat(sourceDir,'tag.txt');
    mlFilename = strcat(sourceDir,'label.txt');

    tmpSImg = load(m1Filename); 
    tmpSDoc = load(m2Filename);
    tmpSLabel = load(mlFilename);
    
    n1Filename = strcat(targetDir,'image.txt');
    n2Filename = strcat(targetDir,'tag.txt');
    nlFilename = strcat(targetDir,'label.txt');

    tmpImg = load(n1Filename);
    tmpDoc = load(n2Filename);
    tmpLabel = load(nlFilename);
    
    %down sampling
    count = histc(tmpSLabel,[1,2]);
    sImg =[]; sDoc = []; sLabel=[];
    for i=1:2
        idx = find(tmpSLabel == i);
        if(count(i) > 3000)
            idx = idx(1:3000,:);
        end
        sImg = [sImg ; tmpSImg(idx,:)];
        sDoc = [sDoc; tmpSDoc(idx,:)];
        sLabel = [sLabel; tmpSLabel(idx,:)];        
    end
    count = histc(tmpLabel,[1,2]);
    tImg =[]; tDoc = []; tLabel=[];
    for i=1:2
        idx = find(tmpLabel == i);
        if(count(i) > 1000)
            idx = idx(1:1000,:);
        end
        tImg = [tImg ; tmpImg(idx,:)];
        tDoc = [tDoc; tmpDoc(idx,:)];
        tLabel = [tLabel; tmpLabel(idx,:)];        
    end
    
    for knn =50
        for sigma = 100:1:100
        	for nc = 200:100:200
                for lambda =1000:100:1000
                    for gamm =100:100:100
                        for sparsity = 0.1
                        
                        dictDirectory = strcat('D:\Dropbox\Dropbox\Project\MSRA\ImgDocResult\Temporary\Dict\car_vs_dog\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamm),'\');
                        
                        ratio = 0.1;
                        rep = 10;
                        N = 100;
%                         
                        [accOur, accB1, accB2, accB3, accB4, accB5] = ExpSCTransferImgDoc(sImg,sDoc,sLabel,tImg,tDoc,tLabel, dictDirectory, ratio, N, rep, sparsity);
                        disp([accB1,accB2,accB3,accB4,accB5, accOur]);
% % % % 
%                         [accB1, accB2] = TestSparsecodingImgDoc(sImg,sDoc,sLabel,tImg,tDoc,tLabel,dictDirectory, ratio, N, rep, sparsity);
%                         result = [result; accB2-accB1];
%                         TestCluster(clusterDirectory,timeslot);
%                         figure()
%                         TestClustering(targetDir,dictDirectory,timeslot,sparsity);
                        end
                    end
                end
            end
        end
    end
end