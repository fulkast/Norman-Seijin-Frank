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
X = double(my_im2col(I, neib));  
M = double(my_im2col(mask, neib));
    
% Construct your dictionary
% If you load your own dictionary U calculated offline you don't have to 
% add anything here
U = buildDictionary(size(X,1));  % TO BE FILLED 
    
% Do the sparse coding with modified Matching Pursuit
Z = sparseCoding(U, X, M, sigma, rc_min);

% Reconstruct the missing parts by using the sparse coding.
n = size(X,2);
for nn = 1:n
    z = Z(:,nn);
    m = M(:,nn);
    u = U(m==0,:);
    x = u*z;
    X(m==0,nn) = x;
end


% Reshape to get I_rec.
I_rec = my_col2im(X, neib, size(I));
