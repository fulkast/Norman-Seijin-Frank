function X = my_im2col(I, patch)

% Provides the functionality of im2col function of the image processing
% toolbox.
%
% INPUT
% I: image
% patch: The size of the square patches extracted
% 
% OUTPUT
% X: (d x n) observations matrix. Obviously d=patch*patch and n is the 
%            number of patches extracted
%
% SPEED
% my_im2col should be about as fast as im2col.

if ndims(I) ~= 2
    error('Only grayscale images supported for now')
end

n = (size(I)'./patch(:));
if ~all(~mod(n, 1))
    warning('Image not divisible by patch size. Boundary pixels will be lost');
end

if numel(patch)==1
    patch = [patch, patch];
end

n  = floor(n);
sx = [patch(1)*patch(2), n(1)*n(2)];

% Preallocate output X (this is way faster than using zeros(sx(1), sx(2)))
X(sx(1), sx(2)) = cast(0, 'like', I);

% % Idea: copy line-by-line, hoping that this goes faster.
% % Block-copy (below) is faster.
% 
% % Rounded image height
% h = n(1)*patch(1);
% for i=1:h
%     for j=1:n(2)
%         % Patch index 
%         p = floor((i-1)/patch(1))*n(2)+j;
%         % Slice begin and end
%         sb = (mod(i-1,patch(1)))*patch(2)+1;
%         se = sb + patch(2)-1;
%         X(sb:se, p) = I(i, (j-1)*patch(2)+1:j*patch(2))';
%     end
% end

% Copy patches to columns of X.
for i=1:n(1)
    for j=1:n(2)
        imin = (i-1)*patch(1)+1;
        imax = min(i*patch(1), size(I,1));
        jmin = (j-1)*patch(2)+1;
        jmax = min(j*patch(2), size(I,2));
        s = min(sx(1), (imax-imin+1)*(jmax-jmin+1));
        X(:,(i-1)*n(2)+j) = reshape(I(imin:imax, jmin:jmax)', s, 1);
    end
end