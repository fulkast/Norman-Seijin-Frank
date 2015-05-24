function a = lassign(X, mu0, Sigma0, mu1, Sigma1)
a = (p(X, mu0, Sigma0) < p(X, mu1, Sigma1)) + 1;
end

% Use the slow version as reference.
function g = p_veryslow(X, mu, Sigma)
d = size(X,1);
g = zeros(size(X,2),1);
coeff = 1/(2*pi)^(d/2)/det(Sigma)^(0.5);
sinv = inv(Sigma);
for i=1:size(X,2)
    D = X(:,i)-mu;
    g(i) = coeff * exp(-0.5*D'*sinv*D);
end
end

% Get rid of the exp per loop and simplify...
function g = p_slow(X, mu, Sigma)
d = size(X,1);
g = zeros(size(X,2),1);
coeff = log(1/det(Sigma)^(0.5));
sinv = inv(Sigma);
for i=1:size(X,2)
    D = X(:,i)-mu;
    g(i) = coeff - 0.5*D'*sinv*D;
end
end

function g = p(X, mu, Sigma)
D = X-repmat(mu,1,size(X,2));
% Slow version, because we calculate a complete nxn matrix and then take
% just the diagonal entries.
% g = log(1/det(Sigma)^0.5) - 0.5*diag(D'/Sigma*D)';

% This is faster
g = log(1/det(Sigma)^0.5) - 0.5*sum(D'/Sigma.*D',2)';
end