function I_rec = inPaintingPyramid(I, mask)
% for time being, we assume image of size 512 512

% Perform the actual inpainting of the image 

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image 
starting=cputime;

% Parameters
ovlp=8; % even, less or equal to half the size of neib
rc_min = 0.01; % rc_min: minimal residual correlation before stopping
neib = 16; % neib: The patch sizes used in the decomposition of the image
sigma = 0.01; % sigma: residual error stopping criterion, normalized by signal norm

shift=[0 0];
I_rec=zeros(size(I));



U = buildDictionary(neib*neib);  % TO BE FILLED 

for l=0:5
i=5-l;
idxmask=find(mask==0);
IC=I;
maskC=mask;
IC(idxmask)=I_rec(idxmask);
maskC(idxmask)=exp(l^2-26)-exp(-26);
maskC=compress(maskC,2^i);
IC=compress(IC,2^i);
n1=size(IC,1);
n2=size(IC,2);


X = my_im2col(IC, neib,ovlp,shift);  
M = my_im2col(maskC, neib,ovlp,shift);  

Z = sparseCodingPyramid(U, X, M, sigma, rc_min);

X_rec=U*Z;
X_rec=X.*M+(1-M).*X_rec;

I_recC=my_col2im(X_rec,neib,[n1,n2],ovlp,shift);

I_rec=expand(I_recC,2^i);

figure
imshow(I_rec);
fprintf('Took this long %.2f \n',cputime-starting)
end
end
% You need to do the image reconstruction using the known image information
% and for the missing pixels use the reconstruction from the sparse coding.
% The mask will help you to distinguish between these two parts.

% TO BE FILLED
