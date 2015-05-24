function I_comp = CompressByChannel(I)

% Result on submission system:
% Settings:
%   I_comp.colorSpace = 'native';
%   dsf = [2., 4., 4.]; % (bug: scale channels differently...)
%   d = [10, 10, 10];
%   I_comp.discretize = true;
%   I_comp.range.target = [0, 255; 0 255; 0 255];
% Results
%   Run-time:             1.444166s
%   Mean squared error:   0.012878
%   Compression rate:     0.039094
%
% Settings:
%   I_comp.colorSpace = 'ycbcr';
%   dsf = [2., 4., 4.];
%   d = [10, 10, 10];
%   I_comp.discretize = true;
%   I_comp.range.target = [0, 255; 0 255; 0 255];
% Results
%   Run-time:             1.554225s
%   Mean squared error:   0.17115   % => ycbcr doesn't pay off...
%   Compression rate:     0.021966
%
% Settings:
%   I_comp.colorSpace = 'native';
%   dsf = [1., 2., 2.]; % (bug: scale channels differently...)
%   d = [10, 10, 10];
%   I_comp.discretize = true;
%   I_comp.range.target = [0, 255; 0 255; 0 255];
% Results
%   Run-time:             3.993983s
%   Mean squared error:   0.0058825
%   Compression rate:     0.087612

%% Settings
%%%%%%%%%%%%%

% Down-sampling.
dsf = [2., 2., 2.];

% Block size.
d = [10, 10, 10];

% Enable discretization.
I_comp.discretize = true;
I_comp.range.target = [0, 255; 0 255; 0 255];

% Image conversion:
% 'native': rgb/gray
% 'ycbcr'
I_comp.colorSpace = 'native';
if size(I,3)==1
    % Force native color space if grayscale.
    I_comp.colorSpace = 'native';
end
%% Calculation
%%%%%%%%%%%%%%%%

switch I_comp.colorSpace
    case 'ycbcr'
        I = rgb2ycbcr(I);
        %I = ycbcr2rgb(I);
end

% Remember original size of image.
[I_comp.orig.h, I_comp.orig.w, I_comp.c] = size(I);
c = I_comp.c;
dsf = dsf(1:c);
d = d(1:c);

% Downsampling
I_comp.dsf = 1.0;
IC = {};
for i=1:c
    IC{i} = imresize(I(:,:,i), 1/dsf(i));
    [I_comp.h(i), I_comp.w(i)] = size(IC{i});
    I_comp.dsf(i) = dsf(i);
end

X = extract(IC,d);

for i=1:c
    [mu, lambda, U] = pca(X{i});
    lambda = diag(lambda);
    k = find_knee(lambda, 5e-2);
    
    if false
        % Code to choose a certain mode and make it visible.
        k = 1;
        U = U(:,k);
        lambda = lambda(k);
    else
        U = U(:,1:k);
    	lambda = lambda(1:k);
    end
    
    Ix = U'*X{i};
    if I_comp.discretize
        mini = min(Ix(:));
        maxi = max(Ix(:));
        minr = I_comp.range.target(i,1);
        maxr = I_comp.range.target(i,2);
        I_comp.range.orig(i,:) = [mini, maxi];
        mm = (maxr-minr)/(maxi-mini);
        qq = minr - mm*mini;
        Ix = uint8(mm*Ix + qq);
    end
    
    I_comp.I{i} = Ix;
    I_comp.lambda{i} = lambda;
    I_comp.U{i} = U;
    I_comp.d(i) = d(i);
end

function XC = extract(IC, d)
assert(all(size(d)==size(IC)), 'Channel mismatch');

% Loop over all channels
c = size(IC,2);
XC = {};

for i=1:c
    % Expand image such that patches fit entirely
    I = IC{i};
    [w,h] = size(I);
    ex_ny = d(i) - mod(h, d(i));
    ex_nx = d(i) - mod(w, d(i));
    
    if ex_ny<d(i)
        I = [ I, repmat(I(:,end), 1, ex_ny)];
    end
    if ex_nx<d(i)
        I = [ I; repmat(I(end,:), ex_nx, 1)];
    end
    
    % Figure out (new) sizes
    [n, m] = size(I);
    assert(mod(n,d(i))==0);
    assert(mod(m,d(i))==0);
    
    fun = @(block_struct) reshape(block_struct.data, 1, size(block_struct.data, 1)*size(block_struct.data, 2));
    XX = blockproc(I, [d(i), d(i)], fun);
    XC{i} = reshape(permute(XX, [2,1]), d(i)*d(i), n/d(i)*m/d(i));
end

function [mu, lambda, U] = pca(X)
% For one channel only

[d,n] = size(X);
mu = mean(X,2);
M  = repmat(mu, 1, size(X,2));
X  = X-M;
sigma = 1/n*X*X';
[U, lambda] = eig(sigma);


function k = find_knee(lambda, error)
lambda = flipud(lambda);
cs = cumsum(lambda);
k = find(cs>error);
if ~isempty(k)
    k = length(lambda)-k(1)+1;
else
    k = 1;
end
k = min(max(k, 5), length(lambda));
