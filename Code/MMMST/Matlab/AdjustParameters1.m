function result = AdjustParameters1()
    i=1;
    for knn =50
        for sigma = 5:1:5
          for timeslot = 23
            result = [];
        	for nc = 500:450:500
                for lambda =500:100:500
                    for gamm =100:100:100
                        for sparsity = 0.1
                            for timeslotD = 6
                        sourceDir = 'D:\Dropbox\Dropbox\Project\MSRA\Result\Features\Beijing\';
                        targetDir = 'D:\Dropbox\Dropbox\Project\MSRA\Result\Features\Shanghai\';
                        
                        clusterDirectory = strcat('D:\Dropbox\Dropbox\Project\MSRA\Result\Temporary\Cluster\Labeled\0.015\100\inter-1\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamm),'\');
                        dictDirectory = strcat('D:\Dropbox\Dropbox\Project\MSRA\Result\Temporary\Dict\Labeled\0.015\100\inter-1\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamm),'\');
                        
%                         clusterDirectory = strcat('D:\Result\Temporary\Cluster\0.015\010\inter-1\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamm),'\');
%                         dictDirectory = strcat('D:\Result\Temporary\Dict\0.015\010\inter-1\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamm),'\');

                        ratio = 0.1;
                        rep = 10;
                        N = 100;
% % % %                         
                        [accOur, accB1, accB2, accB3, accB4, accB5] = ExpSCTransfer(sourceDir, targetDir, dictDirectory, timeslot, timeslotD, ratio, N, rep, sparsity);
                        disp([timeslot, accB1,accB2,accB3,accB4,accB5, accOur]);
% % 
%                         [accB1, accB2] = TestSparsecoding(sourceDir, targetDir, dictDirectory, timeslot, timeslotD, ratio, N, rep, sparsity);
%                          result = [result; accB2-accB1];
%                         TestCluster(clusterDirectory,timeslot);
%                         figure()
%                         TestClustering(targetDir,dictDirectory,timeslot,sparsity);
                                end
                            end
                        end
                    end
                end
%                 save(strcat('result/SH/',num2str(timeslot)),'result');
            end
        end
    end
end