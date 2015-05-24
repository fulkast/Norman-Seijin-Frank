function I_comp = CompressAllChannels(I)

% Results on submission system: 
% downsampling: 2.
% block size:   10.
% Run-time:             1.166679
% Mean squared error:   0.0084671
% Compression rate:     0.18311
% 

% Down-sampling
downsampling = 2.;

% Block size.
d = 10;
%d = floor(d/downsampling);

% Remember original size of image.
[I_comp.orig.h, I_comp.orig.w, I_comp.orig.c] = size(I);

% Downsampling
I_comp.downsampling = 1.0;
if downsampling > 1.0
    I = imresize(I, 1/downsampling);
    I_comp.downsampling = downsampling;
end

% Remeber the new image size.
[I_comp.h, I_comp.w, I_comp.c] = size(I);


X = extract(I,d);
c = size(I,3);
k=-1;

for i=1:c
    [mu, lambda, U] = pca(X(:,:,i));
    lambda = diag(lambda);
    
    if k<0
        k = find_knee(lambda, 5e-2);
    end
    U = U(:,1:k);
    lambda = lambda(1:k);
    
    
    
    I_comp.I(:,:,i) = U'*X(:,:,i); 
    I_comp.lambda(:,i) = lambda;
    I_comp.U(:,:,i) = U;
    I_comp.d = d;
end


function X = extract(I, d)

% Expand image such that patches fit entirely
ex_nx = d - mod(size(I,1), d);
ex_ny = d - mod(size(I,2), d);
if ex_ny<d
    I = [ I, repmat(I(:,end, :), 1, ex_ny)];
end
if ex_nx<d
    I = [ I; repmat(I(end,:, :), ex_nx, 1)];
end

% Figure out sizes
[n, m, c] = size(I);
assert(mod(n,d)==0);
assert(mod(m,d)==0);

%fun = @(block_struct) block_struct.data(:,:,[2 1 3]);
fun = @(block_struct) reshape(block_struct.data, 1, size(block_struct.data, 1)*size(block_struct.data, 2), size(block_struct.data, 3));
X = blockproc(I, [d, d], fun);
X = reshape(permute(X, [2,1,3]), d*d, n/d*m/d, c);

function [mu, lambda, U] = pca(X)
% For one channel only
assert(size(X,3)==1);

[d,n] = size(X);
mu = mean(X,2);
M  = repmat(mu, 1, size(X,2));
X  = X-M;
sigma = 1/n*X*X';
[U, lambda] = eig(sigma);


function k = find_knee(lambda, error)
lambda = flipud(lambda);
k = find(cumsum(lambda)>error);
k = length(lambda)-k(1);

