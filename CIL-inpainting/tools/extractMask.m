
fnMask = 'data/bungee1';
fnOut =  [fnMask, '_mask'];
fnExt = '.png';

color = [0, 255, 0];

I = imread([fnMask, fnExt]);
mask = 255*ones(size(I));
m = I(:,:,1)==color(1) & I(:,:,2)==color(2) & I(:,:,3)==color(3);
mask(m) = 0;
mask = mask(:,:,1);

imwrite(mask, [fnOut, fnExt]);