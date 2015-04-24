function [mu , lambda, U]  = PCAanalyse(X)
[r, ~ ,v] = size(X);

mu = zeros(r,v);
lambda = zeros(r,r,v);
U = zeros(r,r,v);

for i = 1:v

mu(:,i) = mean(X(:,:,i),2);


X(:,:,i) = X(:,:,i) - repmat(mu(:,i),[1 size(X(:,:,i),2)]);

varying = cov(X(:,:,i)');

[ U(:,:,i) , lambda(:,:,i) ] = eig(varying);


end



