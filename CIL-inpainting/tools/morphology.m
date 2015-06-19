% This script assumes that I, mask, I_rec have been loaded.
% (Execute e.g. EvaluateInpainting script)
I = imread('../../data/claudia_512x512.png');
mask = imread('../../data/claudia_512x512_mask.png');
I_rec = imread('../output/claudia_512x512_mask_ref.png');


% Ensure images to be double.
I_rec = im2double(I_rec);
mask = im2double(mask);
I = im2double(I);

% Apply mask to I.
I_mask = I;
I_mask(mask==0) = 0;

% I_mask (and mask) is the input for inpainting
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) Dilate I_mask to get a first guess about the missing pixels.
% Radius defines the shape of the structuring elements.
% radius =3;
% SE = strel('disk', radius);
% X = imdilate(I_mask,SE);
% X(mask~=0) = 0;
filtersize = 2;
H = fspecial('disk',filtersize);
X = imfilter(I_mask, H,'replicate');

if false
    % Try to correctly scale the blur (wihtout artifacts close to corners)
    H = fspecial('average',size(H,1));
    C = imfilter(1-mask, H,'replicate');
    X = X+C;
    %minx = min(X(mask==0));
    %maxx = max(X(mask==0));
    %X = 1/(maxx-minx)*X - minx/(maxx-minx);
    X(mask~=0) = 0;
else
    X = X+0.5;
    X(X>1)=1;
    X = (X - 0.5)*2;
    X(mask~=0) = 0;
end

if 0
    TTT=I_mask;
    TTT(mask==0) = X(mask==0);
    figure; imshow(TTT)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2) Blur the mask that will be used to get a smooth transient...
smoothing = 2;
H = fspecial('disk',smoothing);
Y = imfilter(mask,H,'replicate');
Y = Y+0.5;
Y(Y>1)=1;
Y = (Y - 0.5)*2;
if false
    Y(mask~=0)=1;
    figure(1), imshow(Y);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3) Mix first guess with reconstructed part.
if false
    figure(1), imshow(X)
    figure(2), imshow(Y)
end
R = I_rec;
%R(mask~=0) = 0;

R = X.*(Y) + (1-Y).*R;
toc
if false
    II = I_mask;
    II(mask==0) = X(mask==0);
    figure(1); imshow(II)
    figure(2); imshow(R)
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now do the same, but only where confidence on the result is low.

% 4) The harder the edges along the mask borders
M = conv2(double(mask),[1,1,1;1,-8,1;1,1,1],'same');
M(M<0) = 0;
M(M>0) = 1;
Z = abs(conv2(double(I_rec),[1,1,1;1,-8,1;1,1,1],'same'));
if false
    figure; imshow(M)
    return
end
minZ = min(Z(:)); maxZ = max(Z(:));
Z = 1/(maxZ-minZ)*Z - minZ/(maxZ-minZ);
Z = M.*Z;
threshold = 0.05;
Z(Z<threshold) = 0;
Z(Z>=threshold) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5) Blur the mask that will be used to get a smooth transient...
smoothing = 2;
H = fspecial('disk',smoothing);
Y = imfilter(Z,H,'replicate');
Y = Y+0.5;
Y(Y>1)=1;
Y = (Y - 0.5)*2;
Y(mask~=0)=1;

if true
    figure(1); imshow(Z);
    figure(2); imshow(Y);
    figure(3); imshow(I_rec);
    return
end


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




















% Shift
%s = 0.5;
%xdata=[1 size(I,2)];
%ydata=[1 size(I,1)];
%T=maketform('affine',[1 0 0; 0 1 0; s s 1]);
%C=imtransform(C,T,'bilinear','XData',xdata,'YData',ydata);
