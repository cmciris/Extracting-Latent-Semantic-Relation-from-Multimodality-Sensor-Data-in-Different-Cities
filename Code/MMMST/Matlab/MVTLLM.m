% **********************Baseline: Multi-view Transfer Learning with a Large Margin Approach*********************
% Created by Ying Wei
% 2014/01/30
% Input: X: Source domain input
% 	     Z: Target domain input
% 		 noDim: number of the dimensions of the first view
% 		 labelX: label vector for the X
% Output:
function [accuracy, predLabels] = MVTLLM(X, Z, noDim, labelX, groundTruth)
	options = statset('Display','final');
	gm1 = fitgmdist(X,2,'Options',options);
	gm2 = fitgmdist(Z,2,'Options',options);

	% estimate the ratio betai according to gm1 and gm2
	% 

	param.gamma1 = 0.1;
	param.gamma2 = 0.1;
	param.C1 = 0.1;
	param.C2 = 0.1;
	param.C = 0.1;
	param.epsilo = 0.1;

	X1 = X(:,1:noDim);
	X2 = X(:,(noDim + 1):end);
	Z1 = Z(:,1:noDim);
	Z2 = Z(:,(noDim + 1):end);

	[nX,d1]=size(X1);
	[~,d2]=size(X2);
	[nZ,~]=size(Z1);

	% start to perform the algorithm part
	% construct the Xb, X1b, X2b
	Xb = [X1 -X2]';
	Zb = [Z1 -Z2]';
	X1b = [X1 zeros(nX,d2)]';
	X2b = [zeros(nX,d1) X2]';
	% construct the matrix H(B,C)
	H1 = diag([param.gamma1 * ones(d1,1); param.gamma2/2 * ones(d2,1)]);
	H2 = (beta .* Xb) * Xb' + Zb * Zb';
	H = H1 + param.C * H2;
	% start to iterate
	t = 0; et = 1e+12;
	a = [];
	b = [];
	w = [];
	w = [w rand(d1+d2,1)];
	while(et>param.epsilo)
		t = t+1;
		% calculate the subgradient of Remp(w)
		I1p = labelX .* (w(:,1)' * X1b)';
		I1p(I1p < 1,:) = 1;
		I1p(I1p >=1,:) = 0;
		I2p = labelX .* (w(:,1)' * X2b)';
		I2p(I2p < 1,:) = 1;
		I2p(I2p >=1,:) = 0;
		at = -param.C1 * X1b * (beta .* I1p .* labelX) - param.C2 * X2b * (beta .* I2p .*labelX);
		a = [a at];
		% calculate the Remp(w)
		Remp = param.C1 * sum(beta .* max([zeros(nX,1) ones(nX,1)-labelX .* (w(:,1)' * X1b)'], [], 2))...
			 + param.C2 * sum(beta .* max([zeros(nX,1) ones(nX,1)-labelX .* (w(:,1)' * X2b)'], [], 2));
		bt = Remp - dot(w(:,1),at);
		b = [b bt];
		% derive the optimization problem 
		[wbest,hist]= StochSubg(H,a,b,w(:,1),100);
		% calculate the error
		% calculate the J(wi)
		minJw = [];
		for i=1:t
			theta1 = ones(nX,1) - labelX .* (w(:,i)' * X1b)';
			theta2 = ones(nX,1) - labelX .* (w(:,i)' * X2b)';
			Jw = w(:,i)' * H * w(:,i) + param.C1 * beta .* theta1 + param.C2 * beta .* theta2;
			minJw = [minJw; Jw];
		end
		minJw = min(minJw);
		if(minJw - hist{1} <= param.epsilo)
			break;
		end
	end
	% classification assignment 
	predLabels = -1 * zeros(nZ,1);
	ind = 0.5 * w' * Z' > 0 + 0;
	predLabels(ind', :) = 1;
	accuracy = sum((predLabels == groundTruth) + 0) ./ numel((predLabels == groundTruth) + 0);
end

function [w,whist] = StochSubg(Hbar,abar,bbar,w0,MAX_ITERS)
%********************************************************************
% stochastic subgradient method for expected piecewise linear problem
% uses square summable, but nonsummable step size rule, alpha_k = a/k
%********************************************************************
f = [+Inf]; fbest = [+Inf]; 
iter = 1;
w = w0;
whist = [w];

while iter < MAX_ITERS 
  	if( rem(iter,10) == 0 ), fprintf(1,'iter: %d\n',iter), end

    % find the index corresponding to the maximum value of Rcp
  	tmp = w' * abar + bbar
  	[~,index] = max(tmp,[],2);
  	g = 2 * H * w + abar(:,index);
  	fval = w' * H * w + w' * abar(:,index) + bbar(:,index);

  	% step size selection
 	alpha = 1/iter;

  	% objective values
  	f(end+1) = fval;
  	fbest(end+1) = min( fval, fbest(end) );

  	% subgradient update
  	w = w - alpha*g; iter = iter + 1; whist = [whist, w];
end

% collect history information
hist{1} = fbest; hist{2} = whist;
end
