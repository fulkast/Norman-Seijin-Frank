function I_rec = my_col2im(X, patch, im_size)

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

% It will only work for orthogonal matrices and exact patches coverage

n1=im_size(1);
n2=im_size(2);
d=patch-2;
D=patch;
n1new=(floor(n1/d)+1)*d;
n2new=(floor(n2/d)+1)*d;
Inew=zeros(2+n1new,2+n2new);

%X=zeros(size(I,3)*d*d,n2new*n1new/d^2);
for i=1:n1new/d
    for j=1:n2new/d
        IM=reshape(X(:,(i-1)*n2new/d+j),[D,D,size(X,1)/D^2]);
        Inew((2+1+(i-1)*d):(i*d),(2+1+(j-1)*d):(j*d),:)=IM(3:D-2,3:D-2,:);
        
        Inew((1+(i-1)*d):(2+(i-1)*d),(2+1+(j-1)*d):(j*d),:)=Inew((1+(i-1)*d):(2+(i-1)*d),(2+1+(j-1)*d):(j*d),:)+0.5*IM(1:2,3:D-2);
        Inew((1+i*d):(2+i*d),(2+1+(j-1)*d):(j*d),:)=Inew((1+i*d):(2+i*d),(2+1+(j-1)*d):(j*d),:)+0.5*IM(D-1:D,3:D-2);
        Inew((2+1+(i-1)*d):(i*d),(1+(j-1)*d):(2+(j-1)*d),:)=Inew((2+1+(i-1)*d):(i*d),(1+(j-1)*d):(2+(j-1)*d),:)+0.5*IM(3:D-2,1:2);
        Inew((2+1+(i-1)*d):(i*d),(1+j*d):(2+j*d),:)=Inew((2+1+(i-1)*d):(i*d),(1+j*d):(2+j*d),:)+0.5*IM(3:D-2,D-1:D);
        
        Inew((1+(i-1)*d):(2+(i-1)*d),(1+(j-1)*d):(2+(j-1)*d),:)=Inew((1+(i-1)*d):(2+(i-1)*d),(1+(j-1)*d):(2+(j-1)*d),:)+0.25*IM(1:2,1:2);
        Inew((1+i*d):(2+i*d),(1+(j-1)*d):(2+(j-1)*d),:)=Inew((1+i*d):(2+i*d),(1+(j-1)*d):(2+(j-1)*d),:)+0.25*IM(D-1:D,1:2);
        Inew((1+(i-1)*d):(2+(i-1)*d),(1+(j)*d):(2+(j)*d),:)=Inew((1+(i-1)*d):(2+(i-1)*d),(1+(j)*d):(2+(j)*d),:)+0.25*IM(1:2,D-1:D);
        Inew((1+i*d):(2+i*d),(1+j*d):(2+j*d),:)=Inew((1+i*d):(2+i*d),(1+j*d):(2+j*d),:)+0.25*IM(D-1:D,D-1:D);
    end
end
Inew=2*Inew;
Inew(3:n1new,3:n2new)=Inew(3:n1new,3:n2new)/2;

Inew(1:2,1:2)=Inew(1:2,1:2)*2;
Inew(1:2,n2new-1:n2new)=Inew(1:2,n2new-1:n2new)*2;
Inew(n1new-1:n1new,1:2)=Inew(n1new-1:n1new,1:2)*2;
Inew(n1new-1:n1new,n2new-1:n2new)=Inew(n1new-1:n1new,n2new-1:n2new)*2;

I_rec(:,:,:)=Inew(2:n1+1,2:n2+1,:);