function [accuracy, predLabel] = MDT(X, Z, noDim, labelX, groundTruth)
    idx = labelX == 1;
    labelX(idx,:)=-1;
    labelX(~idx,:)=1;
    idx = groundTruth == 1;
    groundTruth(idx,:)=-1;
    groundTruth(~idx,:)=1;
    idx = sum([X;Z]) == 0;
    
    X(:,idx) = [];
    Z(:,idx) = [];
    noDim = noDim - sum(idx(1:noDim));
    
    ratio = 0.05;
    nR = floor(size(Z,1) * ratio / 2);
    predLabel = [];
    index = [];
    allIndex = [1:size(groundTruth,1)]';
    while(size(Z,1)>0)
        disp([size(Z,1)]);

	X1 = X(:,1:noDim);
	X2 = X(:,(noDim + 1):end);
	Z1 = Z(:,1:noDim);
	Z2 = Z(:,(noDim + 1):end);

	[nX,d1]=size(X1);
	[~,d2]=size(X2);
	[nZ,~]=size(Z1);

	%calculate the Mw
	Mw = X1' * (labelX * labelX') * X2;
	Qw = [zeros(d1,d1) Mw ; Mw' zeros(d2,d2)];

	Mx = cov(X1,1);
	Mz = cov(X2,1);
	P = [Mx zeros(d1,d2); zeros(d2,d1) Mz];

	S1 = [X1;Z1];
	S2 = [X2;Z2];
	L = [ones(nX,nX)/(nX^2) -ones(nX,nZ)/(nX*nZ); -ones(nZ,nX)/(nX*nZ) ones(nZ,nZ)/(nZ^2)];
	Qd = [S1'*L*S1 zeros(d1,d2); zeros(d2,d1) S2'*L*S2];

	Qc = [S1'*S1 -S1'*S2; -S2'*S1 S2'*S2];

	Q = Qw - 16*Qd - 256*Qc;

    P = P + diag(rand(d1+d2,1)*0.00000001); %in case that P is singular
	[V, D] = eigs(Q,P);

	wx = V(1:d1,1)*1000;
	wz = V(d1+1:end,1);

	bx = (mean(X1(labelX == -1,:)) * wx + mean(X1(labelX == 1,:)) * wx) / 2;
	bz = (mean(X2(labelX == -1,:)) * wz + mean(X2(labelX == 1,:)) * wz) / 2;
	b = bx + bz;

	F = Z1 * wx; %+ Z2 * wz - b; 

    [~,n]=sort(F);
    [~,m]=sort(F,'descend');
    if(nR * 2 >= nZ)
        X=[X;Z];
        labelX =[labelX;((F>0)+0)*2-1];
        predLabel = [predLabel;((F>0)+0)*2-1];
        idx = [1:size(Z,1)]';
        Z = [];
        index = [index; allIndex(idx,:)];
        allIndex(idx,:)=[];
    else
    X = [X; Z(n(1:nR),:)];
    X = [X; Z(m(1:nR),:)];
    labelX = [labelX; -((F(n(1:nR))<=0)+0)*2+1];
    labelX = [labelX; ((F(m(1:nR))>=0)+0)*2-1];
    idx = [n(1:nR);m(1:nR)];
    Z(idx,:)=[];
    predLabel = [predLabel;-((F(n(1:nR))<=0)+0)*2+1];
    predLabel = [predLabel;((F(m(1:nR))>=0)+0)*2-1];
    index = [index; allIndex(idx,:)];
	allIndex(idx,:)=[];
    end
    end
    
    [~,sorted]=sort(index);
    predLabel = predLabel(sorted,:);
    accuracy = sum(predLabel == groundTruth) ./ length(groundTruth);
end