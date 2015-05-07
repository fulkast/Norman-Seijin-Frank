function I_rec = inPainting(I, mask)

% Perform the actual inpainting of the image 

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image 

% Parameters
rc_min = 0.01; % rc_min: minimal residual correlation before stopping
neib = 16; % neib: The patch sizes used in the decomposition of the image
sigma = 0.01; % sigma: residual error stopping criterion, normalized by signal norm

% Get patches of size neib x neib from the image and the mask and
% convert each patch to 1D signal
X = my_im2col(I, neib);  
M = my_im2col(mask, neib);  
n1=size(I,1);
n2=size(I,2);

% Construct your dictionary
% If you load your own dictionary U calculated offline you don't have to 
% add anything here
 U = buildDictionary(neib*neib*size(I,3));  % TO BE FILLED 
    
% Do the sparse coding with modified Matching Pursuit
  Z = sparseCoding(U, X, M, sigma, rc_min);
  X_rec=U*Z;
  idx=find(M~=0);
  X_rec(idx)=X(idx);
  I_rec=my_col2im(X_rec,neib,[n1,n2]);
  

% You need to do the image reconstruction using the known image information
% and for the missing pixels use the reconstruction from the sparse coding.
% The mask will help you to distinguish between these two parts.

% TO BE FILLED
