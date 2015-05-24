function Y = stdize(X)
n = size(X,1);
mu = mean(X,1);
stdev = std(X,0,1);
%s = sqrt(1/(n-1)*sum(X.^2,1)-n*m.^2);
Y = (X-repmat(mu, n, 1)) ./ repmat(stdev, n, 1);
end

