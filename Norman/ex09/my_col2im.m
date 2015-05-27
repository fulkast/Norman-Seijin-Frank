function I_rec = my_col2im(X, patch, im_size)

% Provides the functionality of col2im function of the image processing
% toolbox.
%
% INPUT
% X: (d x n) observations matrix. Obviously d=patch*patch and n is the 
%            number of patches extracted
% patch: The size of the square patches extracted
% im_size: Size of the original image 
%
% OUTPUT
% I_rec: image

if numel(patch)>1
    error('Only squared patches are supported');
end

if size(X,1) ~= patch*patch
    error('Patch size does not match the input X');
end

% Number of patches per dimension.
n = floor(im_size/patch);

% Allocate, but don't initialize.
I_rec(im_size(1), im_size(2)) = cast(0, 'like', X);

for j = 1:size(X,2)
    % Patch indices (zero-based).
    pj = mod(j-1, n(2));
    pi = floor((j-1)/n(2));
    pjmin = pj*patch+1;
    pjmax = pjmin+patch-1;
    pimin = pi*patch+1;
    pimax = pimin+patch-1;
    I_rec(pimin:pimax, pjmin:pjmax) = reshape(X(:,j), patch, patch)';
end