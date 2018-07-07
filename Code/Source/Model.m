%%
%%%The below parameters are specified to retrieve the appropriate%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%coupled dictionaries learned%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% OR clustering result%%%%%%%%%%%%%%%%%%%%%%%%%%

% the number of k nearest numbers during constructing the graph
knn = 50;
% the sigma/lambda/gamma parameter during learning sematically related
% dictionaries
sigma = 5;
lambda = 500;
gamma = 100;
% the timeslot of a source domain from which features and dictionaries
% learnt on it are transferred.
timeslotS = 8;
% number of clusters / dictionary size of semantically related dictionaries
nc = 50;


%%
%%%The below parameters are specified for multi-modal transfer adaboost%%%


%the ratio of labelled data used for training in the target domain
ratio = 0.1;
% the timeslot of the target domain of interest
timeslot = 2;
% the number of repetitions of multi-modal transfer adaboost experiments
rep =  10;
% the number of iterations in the multi-modal transfer adaboost algorithm
N = 100;

%%
%%%The directories from which data/dictionaries/clustering results should be retrieved
sourceDir = 'DemoData\source\';
targetDir = 'DemoData\target\';
                        
clusterDirectory = strcat('ModelData\cluster\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamma),'\');
dictDirectory = strcat('ModelData\dictionary\knn_',num2str(knn),'\sigma_',num2str(sigma),'\nc_',num2str(nc),'\lambda_',num2str(lambda),'\gamma_',num2str(gamma),'\');

[accuracy] = ExpSCTransfer(sourceDir, targetDir, dictDirectory, timeslot, timeslotS, ratio, N, rep);

disp(accuracy)