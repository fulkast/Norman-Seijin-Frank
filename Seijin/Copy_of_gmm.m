function [z,Z, U, loglike] = gmm(data, nClusters, varargin)
% GMM   Gaussian Mixture Model with fixed diagonal covariance
% 
%   [z, U, loglike, Z] = gmm(data, nClusters, varargin)
%   Input:
%       data        Data matrix, columns are observations
%       nClusters   Number of clusters
%
%   Output:
%       z           Cluster label vector
%       U           Matrix, columns are estimated cluster means
%       loglike     log-likelihood of the data.
%       Z           responsibilities.
%
%   Additional parameters through varargin
%       threshold   Stop if change of log-likelihood smaller than threshold,
%                   the default is 1e-2.
%       maxIter     Maximum number of iterations, the default is 100.
%       means       Initialization of the cluster centroids.
%       variance    We fix the variance to variance*I, the default is 1e-3.
%       usinglogsumexp
%                   Perform computations in the log domain, the default is 1.

[foo, foo,U, foob]=k_means(data,nClusters);
% initializations
threshold = 1e-2;
maxIter = 20;
[nDims, nExamples] = size(data);
Z = zeros(nClusters, nExamples);
loglike = -realmax;
change = threshold+1;
iterations = 0;
variance = 1e-3;
usingLogSumExp = 1;

% initialization of the class means
%indices = randperm(nExamples);
%U = data(:,indices(1:nClusters));

% uniform prior
pi = (1/nClusters)*ones(nClusters,1);

% parse varargin
for k=1:2:length(varargin)
    switch lower(varargin{k})
    case 'threshold'
        threshold = varargin{k+1};
    case 'maxiter'
        maxIter = varargin{k+1};
    case 'means'
        U = varargin{k+1};
    case 'variance'
        variance = varargin{k+1};
    case 'usinglogsumexp'
        usingLogSumExp = varargin{k+1};
    otherwise
        error(['Unknown parameter ''' varargin{k} '''.']);
    end
end

% initialize covariance matrices
Sigma = zeros(nDims, nDims, nClusters);
for k=1:nClusters
    Sigma(:,:,k) = diag(variance*ones(nDims,1));
end


% GMM algorithm
while (abs(change)>threshold && iterations < maxIter),
    iterations = iterations+1;
 
    % E-step: estimate responsibilities
    for k=1:nClusters
        Z(k,:) = LogGaussPDF(data, U(:,k), Sigma(:,:,k));
    end
    if (usingLogSumExp)
        Z = Z + repmat(log(pi+eps), [1 nExamples]);
        mass = logsumexp(Z,1);
        logP = Z;
        Z = Z - repmat(mass, nClusters, 1);
        Z = exp(Z);
    else
        Z = exp(Z);
        Z = Z.*repmat(pi, [1 nExamples]);
        mass = sum(Z,1);
        P = Z;
        Z = Z./repmat(mass, nClusters, 1);
    end
    
    % M-step: estimate parameters
    U = data*Z';
    U = U./repmat(sum(Z,2)'+eps, nDims, 1);
    pi = sum(Z,2);
    pi = pi./nExamples;
    
    
    Sigma=zeros(nClusters,nDims, nDims);
    
    for n=1:nExamples
        Sigma=Sigma+repmat(Z(:,n),1,nDims,nDims).*permute(repmat(data(:,n)*data(:,n)',1,1,nClusters),[3,1,2]);
    end
    Sigma=Sigma./repmat((sum(Z,2)+eps),1,nDims,nDims);
    Sigma=permute(Sigma,[2,3,1]);
    
    % estimate change of log-likelihood
    loglike_old = loglike;
    if (usingLogSumExp)
        loglike = sum(logsumexp(logP,1));
    else
        loglike = sum(log(sum(P,1)));
    end
    change = loglike - loglike_old;
    
       disp(sprintf('Iterations = %d  Loglikelihood = %0.5g  Change = %0.5g.', ...
                   iterations, loglike, change));

end

% convert assignments to vector representation
[foo, z] = max(Z, [], 1);


%----------------------------------------------------------------------------
function P = LogGaussPDF(data, u, Sigma)
eps=0.001;

Sigma =Sigma +diag(eps*ones(size(Sigma,1),1));

[nDims, nExamples] = size(data);

P = data-repmat(u,1,nExamples);
P = sum(P.*(inv(Sigma)*P), 1);
P = -0.5*P - log(sqrt((2*pi)^nDims * (abs(det(Sigma))+realmin)));


%----------------------------------------------------------------------------
function s = logsumexp(a, dim)
% Returns log(sum(exp(a),dim)) while avoiding numerical underflow.
% Default is dim = 1 (columns).
% logsumexp(a, 2) will sum across rows instead of columns.
% Unlike matlab's "sum", it will not switch the summing direction
% if you provide a row vector.

% Written by Tom Minka
% (c) Microsoft Corporation. All rights reserved.

if nargin < 2
  dim = 1;
end

% subtract the largest in each column
[y, i] = max(a,[],dim);
dims = ones(1,ndims(a));
dims(dim) = size(a,dim);
a = a - repmat(y, dims);
s = y + log(sum(exp(a),dim));
i = find(~isfinite(y));
if ~isempty(i)
  s(i) = y(i);
end
