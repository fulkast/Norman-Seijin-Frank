if false
  fnImg = 'data/claudia_512x512.png';
  fnMask = 'data/claudia_512x512_mask.png';
  fnOut = 'claudia_512x512_mask_g.png';
elseif false
  fnImg = 'data/spiral_512x512.png';
  fnMask = 'data/spiral_512x512_mask.png';
  fnOut = 'spiral_512x512_mask_g.png';
elseif true
  fnImg = 'data/TomAndJerry_512x512.png';
  fnMask = 'data/TomAndJerry_512x512_mask.png';
  fnOut = 'TomAndJerry_512x512_mask_g.png';
end

I = imread(fnImg);
mask = imread(fnMask);

Jr = I;
Jg = I;
Jb = I;
Jr(mask==0) = 0;
Jg(mask==0) = 255;
Jb(mask==0) = 0;

J(:,:,1) = Jr;
J(:,:,2) = Jg;
J(:,:,3) = Jb;

imwrite(J, fnOut);