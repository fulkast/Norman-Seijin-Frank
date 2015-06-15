% This script assumes that I, mask, I_rec have been loaded.
% (Execute e.g. EvaluateInpainting script)

% Apply mask to I.
I_mask = I;
I_mask(mask==0) = 0;

% I_mask (and mask) is the input for inpainting

% 1) Dilate I_mask to get a first guess about the missing pixels.
% Radius defines the shape of the structuring elemente.
tic
radius =3;
SE = strel('disk', radius);
X = imdilate(I_mask,SE);
X(mask~=0) = 0;

% 2) Blur the mask that will be used to get a smooth transient...
smoothing = 2;
H = fspecial('disk',smoothing);
Y = imfilter(mask,H,'replicate');
Y = Y+0.5;
Y(Y>1)=1;
Y = (Y - 0.5)*2;
%Y(mask~=0)=1;
%imshow(Y);


% 3) Mix first guess with reconstructed part.
R = I_rec;
R(mask~=0) = 0;

R = X.*(1-Y) + Y.*R;
II = I_mask;
II(mask==0) = R(mask==0);
toc


% Now do the same, but only where confidence on the result is low.

% 4) The harder the edges along the mask borders
M = conv2(double(mask),[1,1,1;1,-8,1;1,1,1],'same');
M(M>0) = 1;
M(M<1) = 0;
Z = abs(conv2(double(I_rec),[1,1,1;1,-8,1;1,1,1],'same'));
minZ = min(Z(:)); maxZ = max(Z(:));
Z = 1/(maxZ-minZ)*Z - minZ/(maxZ-minZ);
Z = M.*Z;
threshold = 0.1;
Z(Z<threshold) = 0;
Z(Z>=threshold) = 1;

% 5) Blur the mask that will be used to get a smooth transient...
smoothing = 2;
H = fspecial('disk',smoothing);
Y = imfilter(Z,H,'replicate');
Y = Y+0.5;
Y(Y>1)=1;
Y = (Y - 0.5)*2;
Y(mask~=0)=1;

% 6) Mix first guess with reconstructed part.
R = I_rec;
%R(mask~=0) = 0;

R = X.*(1-Y) + Y.*R;
II = I_rec;
II(Y~=0) = R(Y~=0);
toc
figure(1); imshow(II)
figure(2); imshow(I_rec)

Error1 = mean(mean(mean(((I - I_rec) ).^2)))
Error2 = mean(mean(mean(((I - II) ).^2)))

% Figure



