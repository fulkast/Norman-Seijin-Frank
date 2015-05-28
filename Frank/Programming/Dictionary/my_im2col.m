function X = my_im2col(I, patch)

% Provides the functionality of im2col function of the image processing
% toolbox.
%
% INPUT
% I: image
% patch: The size of the square patches extracted
% 
% OUTPUT
% X: (d x n) observations matrix. Obviously d=patch*patch and n is the 
%            number of patches extracted

% You can write a for-loop to extract the patches one by one and then 
% transform each patch to an 1D signal sequentially
% creating matrix X

n1=size(I,1);
n2=size(I,2);
d=patch-2;
D=patch;
n1new=(floor(size(I,1)/d)+1)*d;
n2new=(floor(size(I,2)/d)+1)*d;
Inew=zeros(n1new+2,n2new+2,size(I,3));
Inew(2:n1+1,2:n2+1,:)=I(:,:,:);
X=zeros(size(I,3)*D*D,n2new*n1new/d^2);
for i=1:n1new/d
    for j=1:n2new/d
        X(:,(i-1)*n2new/d+j)=reshape(Inew((1+(i-1)*d):(2+i*d),(1+(j-1)*d):(2+j*d),:),[size(I,3)*D*D,1]);
    end
end
end