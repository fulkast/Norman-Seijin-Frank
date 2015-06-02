function I_rec = my_col2im(X, patch, im_size,ovlp,shift)
shift(1)=mod(shift(1),patch-ovlp);
shift(2)=mod(shift(2),patch-ovlp);
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
d=patch-ovlp;
D=patch;
n1new=(ceil((n1+shift(1))/d))*d;
n2new=(ceil((n2+shift(2))/d))*d;
Inew=zeros(ovlp+n1new,ovlp+n2new);

%X=zeros(size(I,3)*d*d,n2new*n1new/d^2);
for i=1:n1new/d
    for j=1:n2new/d
        
        IM=reshape(X(:,(i-1)*n2new/d+j),[D,D,size(X,1)/D^2]);
        for k=1:ovlp/2
            IM(:,k)=IM(:,k)*k/ovlp;
            IM(k,:)=IM(k,:)*k/ovlp;
            IM(:,D-k+1)=IM(:,D-k+1)*k/ovlp;
            IM(D-k+1,:)=IM(D-k+1,:)*k/ovlp;
        end
        for k=ovlp/2+1:ovlp
            IM(:,k)=IM(:,k)*(k-1)/ovlp;
            IM(k,:)=IM(k,:)*(k-1)/ovlp;
            IM(:,D-k+1)=IM(:,D-k+1)*(k-1)/ovlp;
            IM(D-k+1,:)=IM(D-k+1,:)*(k-1)/ovlp;
        end        
        
        
        
        Inew((1+(i-1)*d):(ovlp+i*d),(1+(j-1)*d):(ovlp+j*d),:)=Inew((1+(i-1)*d):(ovlp+i*d),(1+(j-1)*d):(ovlp+j*d),:)+IM;

    end
end


for k=ovlp/2+1:ovlp
    Inew(:,k)=Inew(:,k)*ovlp/(k-1);
    Inew(:,ovlp+n2new-k+1)=Inew(:,ovlp+n2new-k+1)*ovlp/(k-1);
    Inew(k,:)=Inew(k,:)*ovlp/(k-1);
    Inew(ovlp+n1new-k+1,:)=Inew(ovlp+n1new-k+1,:)*ovlp/(k-1);
end
I_rec(:,:,:)=Inew(1+ovlp/2+shift(1):shift(1)+n1+ovlp/2,shift(2)+1+ovlp/2:shift(2)+n2+ovlp/2,:);



