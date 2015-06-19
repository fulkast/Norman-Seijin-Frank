function [I_rec, stats] = inPainting(I, mask, settings)
% Perform the actual inpainting of the image

% INPUT
% I: (n x n) masked image
% mask: (n x n) the mask hidding image information
%
% OUTPUT
% I_rec = Reconstructed image

%% Parameters
overlap=8;      % even, less or equal to half the size of neib
rc_min = 0.1;  % rc_min: minimal residual correlation before stopping
neib = 16;      % neib: The patch sizes used in the decomposition of the image
sigma = 0.01;   % sigma: residual error stopping criterion, normalized by signal norm
verbose = true; %
advancedBlend = true;

if nargin == 3
    overlap = settings.overlap;
    rc_min = settings.rc_min;
    neib = settings.neib;
    sigma = settings.sigma;
    verbose = settings.verbose;
    advancedBlend = settings.advancedBlend;
end

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
stats.sparsity = numel(Z(Z==0))/numel(Z);
profiling(end+1:end+2) = {cputime-starttime, 'Sparse coding'};

% Reconstruct the missing pixels using the sparse coding.
X_rec=U*Z;
idx=find(M~=0);
X_rec(idx)=X(idx);

if ~advancedBlend
    I_rec=my_col2im(X_rec,neib,[n1,n2],overlap,shift,true);
else
    I_rec_1=my_col2im(X_rec,neib,[n1,n2],overlap,shift,10);
    I_rec_2=my_col2im(X_rec,neib,[n1,n2],overlap,shift,20);
    I_rec_3=my_col2im(X_rec,neib,[n1,n2],overlap,shift,30);
    I_rec_4=my_col2im(X_rec,neib,[n1,n2],overlap,shift,40);
    % figure(1); imshow(I_rec);
    % figure(2); imshow(I_rec_1);
    % figure(3); imshow(I_rec_2);
    % figure(4); imshow(I_rec_3);
    % figure(5); imshow(I_rec_4);
    [I_rec, statistics] = patch_selector(I_rec_1, I_rec_2, I_rec_3, I_rec_4, mask, neib);
    %[I_rec, statistics] = patch_selector2(I_rec_1, I_rec_2, I_rec_3, I_rec_4, mask, neib);
    if verbose
        fprintf('Average number of patches used for blending %g\n', statistics.avg_nop);
    end
end

runtime = cputime-starttime;
profiling(end+1:end+2) = {runtime , 'Reconstruct'};

stats.profiling = profiling;
stats.runtime = runtime;

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