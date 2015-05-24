function I_comp = CompressByChannelSVD(I)

% Result on submission system:
% 

%% Settings
%%%%%%%%%%%%%

% Block size.
d = [15, 15, 15];

% Tolerated error.
error = 0.2;

% Image conversion:
% 'native': rgb/gray
% 'ycbcr'
I_comp.colorSpace = 'native';
if size(I,3)==1
    % Force native color space if grayscale.
    I_comp.colorSpace = 'native';
end

% Down-sampling.
switch I_comp.colorSpace
    case 'native'
        dsf = [1., 1., 1.];
    case 'ycbcr'
        dsf = [1., 2., 2.];
    otherwise
        assert(false, 'Choose the down sampling for this color space')
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

for i=1:c
    % Approximation of X{i} through SVD.
    [U, S, V] = svd(IC{i});
    sv = diag(S);
    k = find_knee(sv, error);
    if false
        % Code to choose a certain mode.
        k = 1;
        scale = 10*[1, 1, 1]; % To make things better visible.
        V = V(:,k);
        U = U(:,k);
        sv = scale(i)*sv(k); 
    else
        V = V(:,1:k);
        U = U(:,1:k);
        sv = sv(1:k);
    end

    I_comp.sv{i} = sv;
    I_comp.U{i} = U;
    I_comp.V{i} = V;
    I_comp.d(i) = d(i);
end

function k = find_knee(lambda, error)
%% Find the "knee" for the sorted lambdas
% New: with error relative to the integral.
lambda = flipud(lambda);
cs = cumsum(lambda);
integral = cs(end);
k = find(cs/integral>error);
if ~isempty(k)
    k = length(lambda)-k(1)+1;
else
    k = 1;
end
k = min(max(k, 5), length(lambda));
