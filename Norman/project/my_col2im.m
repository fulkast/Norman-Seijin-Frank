function I_rec = my_col2im(X, patch, imSize , overlap, shift)

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

if nargin < 5
    shift = [0, 0];
end
if nargin < 4
    overlap = 0;
end

if numel(patch)==1
    patch = [patch, patch];
end

if numel(overlap)==1
    overlap = [overlap, overlap];
end

% if size(X,1) ~= prod(n)
%     error('Patch size does not match the input X');
% end

if any(shift > 0)
    error('Shifts are not supported yet');
    % Need to adapt for loop.
end

% Clip shift.
shift(1)=mod(shift(1),patch(1)-overlap(1));
shift(2)=mod(shift(2),patch(2)-overlap(1));

% Patch sizes.
d = patch - overlap;
D = patch;

% Number of patches per dimension.
n = ceil([(imSize(1)+shift(1))/d(1), (imSize(2)+shift(2))/d(2)]);

% Number of channels.
if numel(imSize) == 3
    c = imSize(3);
else
    c = 1;
end

% If image has to be shifted.
nNew = n.*d+overlap;

% Allocate, but don't initialize.
I_rec(nNew(1), nNew(2), c) = cast(0, 'like', X);

for j = 1:size(X,2)
    % Patch indices (zero-based).
    pj = mod(j-1, n(2));
    pi = floor((j-1)/n(2));
    pimin = pi*d(1)+1;
    pimax = pimin+d(1)-1;
    pjmin = pj*d(2)+1;
    pjmax = pjmin+d(2)-1;
    if any(overlap > 0)
        pp = reshape(X(:,j), D(1), D(2));
        I_rec(pimin:pimax, pjmin:pjmax) = pp(1:d, 1:d);
    else
        % This is faster.
        I_rec(pimin:pimax, pjmin:pjmax) = reshape(X(:,j), D(1), D(2));
    end
end

% If image has to be shifted.
if any(shift > 0) || any(overlap > 0)
    imin = shift(1)+1+overlap(1)/2;
    imax = shift(1)+imSize(1)+overlap(1)/2;
    jmin = shift(2)+1+overlap(2)/2;
    jmax = shift(2)+imSize(2)+overlap(2)/2;
    I_rec = I_rec(imin:imax,jmin:jmax,:);
end


