% *********************multi-modality sparse coding***********************
% Created by Ying Wei
% 2014/12/29
% Inputs: lambda: parameter for sparse coding, controlling the L2
%         pooling: boolean value to determine whether to perform pooling
%         reconstruction error
%         varargin: D,X for different modalities
% Outputs: A: pooling = 1 unified sparse codes by max pooling
%             pooling = 2 unified sparse codes by integrated sc
%             pooling = 0 sparse codes in a matrix for each modality
function [A] = MMSparseCoding(lambda, pooling, varargin)
addpath spams-matlab
start_spams;
% the number of total modalities 
numModalities = (nargin - 2)/2;
% declare the sparse codes container
Alpha = [];
% start to sparse coding for each modality
if(pooling~=2)
    for i=1:numModalities
        D = varargin{i * 2 - 1};
        X = varargin{i * 2};
        %config the lars algorithm's parameter
        param.lambda = lambda(i);
        param.numThreads = -1;
        param.mode = 2;
    
        tic
        alpha = mexLasso(X',D',param);
        alpha = full(alpha);
        alpha(~any(D,2),:) = alpha(~any(D,2),:) - 1e12;
        t=toc
        fprintf('%f instances processed per second\n',size(X,2)/t);
        Alpha = [Alpha; alpha];
    end
    Alpha(Alpha == -1e+12) = 0;
end
% max pooling on the sparse codes 
if(pooling == 1)
    if(strcmp(varargin{end},'imgdoc'))
        Alpha = FeaturePreprocess(Alpha',2); Alpha = Alpha';
    else
%         Alpha = FeaturePreprocess(Alpha',0); Alpha = Alpha';
    end
    p = size(D,1);
    A=[];
    for i =1:p
        tmp = [];
        for j=1:numModalities
            t = Alpha(p*(j-1)+i,:);
            tmp = [tmp;t];
        end
        A = [A;max(tmp)];
    end
elseif(pooling == 2)
    MX = []; 
    MD = []; 
    for i=1:numModalities
        D = varargin{i * 2 - 1};
        X = varargin{i * 2};
        MX = [MX;X'];
        MD = [MD;D'];
    end
    param.lambda = lambda;
    param.numThreads = -1;
    param.mode = 2;
    
    tic;
    A = mexLasso(MX,MD,param);
    A = full(A);
    t=toc;
    fprintf('%f instances processed per second\n',size(MX,2)/t);
    A = FeaturePreprocess(A',2); A = A';
else
    A = Alpha;
end
end