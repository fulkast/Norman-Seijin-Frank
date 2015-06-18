function showDico
dim=16;
t=load('dictionary.mat');
U=t.U;
l=size(U,2);
Usize=[dim*sqrt(l) dim*sqrt(l)]

figure
imshow(0.5+my_col2im(0.5*(U-repmat((mean(U)),dim^2,1))./repmat(sqrt(var(U)),dim^2,1),dim,Usize,0,[0 0]));

end