function X = FeaturePreprocess(X, type)
    X(isnan(X))=0;
    X(isinf(X))=0;
	[nX,~]=size(X);
	if(type == 0)  %
		X = (X-repmat(mean(X),nX,1)) ./ repmat(std(X),nX,1);
	elseif(type == 1) % scaling
		X = (X-repmat(min(X),nX,1)) ./ (repmat(max(X),nX,1) - repmat(min(X),nX,1));
	else %normalization
		for i =1:nX
			X(i,:) = X(i,:) / norm(X(i,:));
		end
    end
    X(isnan(X)) = 0;
end