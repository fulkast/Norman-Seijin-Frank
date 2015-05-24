function I_rec = DecompressByChannelSVD(I_comp)
w = I_comp.w;
h = I_comp.h;
c = I_comp.c;
d = I_comp.d;

assert(length(w)==c);
assert(length(h)==c);
assert(length(d)==c);

I_rec = zeros(I_comp.orig.h, I_comp.orig.w, c);
for i=1:c
    I = I_comp.U{i}*diag(I_comp.sv{i})*I_comp.V{i}'; 
    I_rec(:,:,i) = imresize(I, [I_comp.orig.h, I_comp.orig.w]);
end

switch I_comp.colorSpace
    case 'ycbcr'
        I_rec = ycbcr2rgb(I_rec);
end