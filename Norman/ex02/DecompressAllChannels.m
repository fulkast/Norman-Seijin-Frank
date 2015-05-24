function I_rec = DecompressAllChannels(I_comp)

w = I_comp.w;
h = I_comp.h;
c = I_comp.c;
d = I_comp.d;

npx = ceil(w/d);
npy = ceil(h/d);
I_rec = zeros(npy*d,npx*d,c);

for i=1:c
    X = I_comp.U(:,:,i)*I_comp.I(:,:,i); 
    assert(size(X,1)==d*d);
    n = size(X,2);
    I_rec(1:d, 1:d, i) = reshape(X(:,1), d, d);
    for j=1:n
        px = mod(j-1, npx);
        py = floor((j-1)/npx);
        xl = 1+px*d;
        yl = 1+py*d;
        I_rec(yl:yl+d-1, xl:xl+d-1, i) = reshape(X(:,j), d, d);
    end
end
I_rec = I_rec(1:h,1:w,:);

% If downsampling: increase image.
if I_comp.downsampling > 1.0
    I_rec = imresize(I_rec, [I_comp.orig.h, I_comp.orig.w]);
end
