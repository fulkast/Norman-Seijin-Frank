function [mu, lambda, U]=PCAnalyse(X)
mu=mean(X');
Xx=X-repmat(mu',1,size(X,2));
Xx=Xx-repmat(mean(Xx),size(X,1),1);
[U,lambda]=eig(cov(Xx'));
%diag(lambda)
end