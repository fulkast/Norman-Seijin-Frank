function [z,Z, U, score] = k_means(data, nClusters, varargin)

% K_MEANS   k-means clustering
% 
%   [z, U, score] = k_means(data, nClusters, varargin)
%   Input:
%       data        Data matrix, columns are observations
%       nClusters   Number of clusters
%
%   Output:
%       z           Cluster label vector
%       U           Matrix, columns are estimated cluster means
%       score       Negative log-likelihood of the estimated solution
%                   given the data.
%
%   Additional parameters through varargin
%       threshold   Stop if change of log-likelihood smaller than threshold,
%                   the default is 1e-4.
%       maxIter     Maximum number of iterations, the default is 20.
%       means       Initialization of the cluster centroids.


% initializations
threshold = 1e-4;
maxIter = 30;
[nDims, nExamples] = size(data);
Z = zeros(nClusters, nExamples);
score = realmax;
change = threshold+1;
iterations = 0;

% initialization of the class means
indices = randperm(nExamples);
U = data(:,indices(1:nClusters));

% parse varargin
for k=1:2:length(varargin)  
    switch lower(varargin{k})
    case 'threshold'
        threshold = varargin{k+1};
    case 'maxiter'
        maxIter = varargin{k+1};
    case 'means'
        U = varargin{k+1};
    otherwise
        error(['Unknown parameter ''' varargin{k} '''.']);
    end
end

% k-Means algorithm
while (change>threshold && iterations < maxIter),
    iterations = iterations+1;
   disp(sprintf('Iterations = %d  Change = %0.5g.', iterations, change));

    % E-step: estimate class indices
    for k=1:nClusters,
        Z(k,:) = sum((data-repmat(U(:,k),1,nExamples)).^2,1);
    end
    % assign to maximum
    [foo,ind] = min(Z, [], 1);
    Z = zeros(size(Z));
    Z(sub2ind(size(Z), ind, 1:nExamples)) = 1;
    
    % M-step: estimate means
    U = data*Z';
    U = U./repmat(sum(Z,2)'+eps,nDims,1);

    % estimate change
    score_old = score;
    score = sum(sum((data - U*Z).^2, 1));
    change = score_old - score;
end

% convert assignments to vector representation
[foo, z] = max(Z, [], 1);
