function X = my_im2col(I, patch, overlap, shift)

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
%
% NOTE
% Added shift and overlap

if nargin < 4
    shift = [0, 0];
end
if nargin < 3
    overlap = 0;
end

if numel(patch)==1
    patch = [patch, patch];
end

if numel(overlap)==1
    overlap = [overlap, overlap];
end

if any(overlap > patch/2)
    error('Overlap cannot be larger than half the patch size');
end

% Clip shift.
shift(1)=mod(shift(1),patch(1)-overlap(1));
shift(2)=mod(shift(2),patch(2)-overlap(1));

% Patch sizes.
d = patch - overlap;
D = patch;

% Number of patches per dimension.
n = ceil([(size(I,1)+shift(1))/d(1), (size(I,2)+shift(2))/d(2)]);

% Size of the output matrix.
sx = [size(I,3)*D(1)*D(2), n(1)*n(2)];

% If image has to be shifted.
if any(shift > 0) || any(overlap > 0)
    nNew = n.*d+overlap;
    Inew = zeros(nNew(1), nNew(2), size(I,3));
    imin = shift(1)+1+overlap(1)/2;
    imax = shift(1)+size(I,1)+overlap(1)/2;
    jmin = shift(2)+1+overlap(2)/2;
    jmax = shift(2)+size(I,2)+overlap(2)/2;
    Inew(imin:imax,jmin:jmax,:)=I(:,:,:);
    I = Inew;
end

% Preallocate output X (this is way faster than using zeros(sx(1), sx(2)))
X(sx(1), sx(2)) = cast(0, class(I));

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
        imin = (i-1)*d(1)+1;
        imax = min(i*d(1)+overlap(1), size(I,1));
        jmin = (j-1)*d(2)+1;
        jmax = min(j*d(2)+overlap(2), size(I,2));
        s = min(sx(1), (imax-imin+1)*(jmax-jmin+1)*size(I,3));
        X(:,(i-1)*n(2)+j) = reshape(I(imin:imax, jmin:jmax, :), s, 1);
    end
end

end