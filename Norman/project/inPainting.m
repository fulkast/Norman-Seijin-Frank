function I_rec = inPainting(I, mask)
% Perform the actual inpainting of the image

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image

%% Parameters
overlap=8;      % even, less or equal to half the size of neib
rc_min = 0.01;  % rc_min: minimal residual correlation before stopping
neib = 16;      % neib: The patch sizes used in the decomposition of the image
sigma = 0.01;   % sigma: residual error stopping criterion, normalized by signal norm
verbose = true; %

%% Go!
shift=[0 0];
n1=size(I,1);
n2=size(I,2);
starttime = cputime;
profiling = {};

% Get dictionary. Load it from dictionary.mat if available.
U = buildDictionary(neib*neib);
profiling(end+1:end+2) = {cputime-starttime, 'Build dictionary'};

% Extract patches from images.
M = my_im2col(mask, neib, overlap, shift);
X = my_im2col(I,    neib, overlap, shift);
profiling(end+1:end+2) = {cputime-starttime, 'Extract patches'};

% Sparse coding for the known (unmasked) parts of the image.
Z = sparseCoding(U, X, M, sigma, rc_min);
profiling(end+1:end+2) = {cputime-starttime, 'Sparse coding'};

% Reconstruct the missing pixels using the sparse coding.
X_rec=U*Z;
idx=find(M~=0);
X_rec(idx)=X(idx);
I_rec=my_col2im(X_rec,neib,[n1,n2],overlap,shift);
profiling(end+1:end+2) = {cputime-starttime, 'Reconstruct'};

if verbose
    fprintf('Profiling:\n');
    fprintf('----------------------------\n');
    tstart = 0;
    for i=1:2:length(profiling)
        tend = profiling{i};
        tdelta = tend - tstart;
        fprintf('   %s: %gs\n', profiling{i+1}, tdelta);
        tstart = tend;
    end
    fprintf('   Total: %gs\n', profiling{end-1});
    fprintf('----------------------------\n');
end

end
