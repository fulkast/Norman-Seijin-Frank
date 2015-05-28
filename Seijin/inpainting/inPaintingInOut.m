function I_rec = inPaintingInOut(I, mask)

% Perform the actual inpainting of the image 

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image 
n1=size(I,1);
n2=size(I,2);
toFill=n1*n2;


% Parameters
rc_min = 0.01; % rc_min: minimal residual correlation before stopping
neib = 16; % neib: The patch sizes used in the decomposition of the image
sigma = 0.01; % sigma: residual error stopping criterion, normalized by signal norm

condition=[0 1 0;1 1 1;0 1 0];
 U = buildDictionary(neib*neib);  % TO BE FILLED 
 X = my_im2col(I, neib);
    
while (toFill>0)

    maskSparse=filterImage(mask,condition);
idx=find(maskSparse~=0);
maskSparse(idx)=1;
toFill=n1*n2-numel(idx);
% Construct your dictionary
% If you load your own dictionary U calculated offline you don't have to 
% add anything here

M = my_im2col(mask, neib);  
MS= my_im2col(maskSparse, neib);  
I_rec=zeros(size(I));
 % Get patches of size neib x neib from the image and the mask and
% convert each patch to 1D signal

  
    % Do the sparse coding with modified Matching Pursuit
      Z = sparseCoding(U, X, M, MS,sigma, rc_min);
      X_rec=U*Z;
      idx=find(M~=0);
      X_rec(idx)=X(idx);
    X=X_rec.*MS;
    mask=my_col2im(MS,neib,[n1,n2]);
    I_rec=my_col2im(X,neib,[n1,n2]);
    imshow(I_rec);
end

  

% You need to do the image reconstruction using the known image information
% and for the missing pixels use the reconstruction from the sparse coding.
% The mask will help you to distinguish between these two parts.

% TO BE FILLED
