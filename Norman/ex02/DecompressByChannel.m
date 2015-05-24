function I_rec = DecompressByChannel(I_comp)

w = I_comp.w;
h = I_comp.h;
c = I_comp.c;
d = I_comp.d;

assert(length(w)==c);
assert(length(h)==c);
assert(length(d)==c);

npx = ceil(w./d);
npy = ceil(h./d);

I_rec = zeros(I_comp.orig.h, I_comp.orig.w, c);

for i=1:c
    I = zeros(npy(i)*d(i),npx(i)*d(i));
    Ic = I_comp.I{i};
    if I_comp.discretize
        mini = I_comp.range.orig(i,1);
        maxi = I_comp.range.orig(i,2);
        minr = I_comp.range.target(i,1);
        maxr = I_comp.range.target(i,2);
        mm = (maxi-mini)/(maxr-minr);
        qq = mini - mm*minr;
        Ic = double(Ic);
        Ic = mm*Ic + qq;
    end
    X = I_comp.U{i}*Ic; 
    assert(size(X,1)==d(i)*d(i));
    n = size(X,2);
    for j=1:n
        px = mod(j-1, npx(i));
        py = floor((j-1)/npx(i));
        xl = 1+px*d(i);
        yl = 1+py*d(i);
        I(yl:yl+d(i)-1, xl:xl+d(i)-1) = reshape(X(:,j), d(i), d(i));
    end
    
    I = I(1:h(i),1:w(i),:);
    I_rec(:,:,i) = imresize(I, [I_comp.orig.h, I_comp.orig.w]);
end

switch I_comp.colorSpace
    case 'ycbcr'
        I_rec = ycbcr2rgb(I_rec);
end