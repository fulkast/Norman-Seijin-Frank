function I_rec = inPainting(I, mask)


% Perform the actual inpainting of the image 

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image 

% Parameters
ovlp=8; % even, less or equal to half the size of neib
rc_min = 0.01; % rc_min: minimal residual correlation before stopping
neib = 16; % neib: The patch sizes used in the decomposition of the image
sigma = 0.0001; % sigma: residual error stopping criterion, normalized by signal norm

condition=[0 1 0;1 1 1;0 1 0];
shift=[0 0];
n1=size(I,1);
n2=size(I,2);

% while (numel(find(mask==0)~=0))
maskSparse=filterImage(mask,condition);
idx=find(maskSparse~=0);
maskSparse(idx)=1;
% Construct your dictionary
% If you load your own dictionary U calculated offline you don't have to 
% add anything here
 U = buildDictionary(neib*neib);  % TO BE FILLED 
 
 
 
 
M = my_im2col(mask, neib,ovlp,shift);  
MS= my_im2col(maskSparse, neib,ovlp,shift);  
I_rec=zeros(size(I));
 % Get patches of size neib x neib from the image and the mask and
% convert each patch to 1D signal
% if (size(I,3)==3)
%     %I=rgb2hsv(I);
%     for i=1:3
%           X = my_im2col(I(:,:,i), neib,ovlp);
%           Z = sparseCoding(U, X, M,MS, sigma, rc_min);
%           X_rec=U*Z;
%           idx=find(M~=0);
%           %X_rec(idx)=X(idx);
%           I_rec(:,:,i)=my_col2im(X_rec,neib,[n1,n2],ovlp);
%     end
%    % I_rec=hsv2rgb(I_rec);
%     
% else

    % Do the sparse coding with modified Matching Pursuit


    X = my_im2col(I, neib,ovlp,shift);  
      [Z,newM] = sparseCoding(U, X, M, MS,sigma, rc_min);
      X_rec=U*Z;
      idx=find(M~=0);
      X_rec(idx)=X(idx);
     mask=my_col2im(newM,neib,[n1,n2],ovlp,shift);
      I_rec=my_col2im(X_rec,neib,[n1,n2],ovlp,shift);
     I=I_rec;

   mask(find(mask~=1))=0;
     shift(1)=shift(1)+(neib-ovlp)/2;
     shift(2)=shift(2)+(neib-ovlp)/2;
      ovlp=ovlp*2;
      neib=neib*2;
%    imshow(mask.*I);
end
%end
  

% You need to do the image reconstruction using the known image information
% and for the missing pixels use the reconstruction from the sparse coding.
% The mask will help you to distinguish between these two parts.

% TO BE FILLED
