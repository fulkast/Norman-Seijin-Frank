function GMMlearning()

foo=load('DATA.mat');
X=foo.X;

 [z,Z, U]=gmm(X, 100);

Y=U*Z;
imshow((my_col2im(Y,16,[512 512],0,[0 0]))  );