function X = my_im2col(I, patch,ovlp,shift)
shift(1)=mod(shift(1),patch-ovlp);
shift(2)=mod(shift(2),patch-ovlp);
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
d=patch-ovlp;%ovlp<=patch/2
D=patch;
n1new=(ceil((size(I,1)+shift(1))/d))*d;
n2new=(ceil((size(I,2)+shift(2))/d))*d;
% Inew=-0.01*ones(n1new+ovlp,n2new+ovlp,size(I,3));
Inew=zeros(n1new+ovlp,n2new+ovlp,size(I,3));
Inew(shift(1)+1+ovlp/2:shift(1)+n1+ovlp/2,shift(2)+1+ovlp/2:shift(2)+n2+ovlp/2,:)=I(:,:,:);
X=zeros(size(I,3)*D*D,n2new*n1new/d^2);
for i=1:n1new/d
    for j=1:n2new/d
        X(:,(i-1)*n2new/d+j)=reshape(Inew((1+(i-1)*d):(ovlp+i*d),(1+(j-1)*d):(ovlp+j*d),:),[size(I,3)*D*D,1]);
    end
end
end