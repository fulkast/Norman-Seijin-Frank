function [I_rec, stats] = patch_selector(I1, I2, I3, I4, mask, patch)
% stats optionally informs caller about some statistics
%   stats.avg_nop:      average number of patches used for blending
%   stats.confidences:  confidence data per patch

CONFIDENCE_THRESHOLD = 1.5;
ENABLE_LOCAL_WEIGHING = true;
G_SIZE = 3;
G_SIGMA = 1;

II = {I1, I2, I3, I4};

if numel(patch) == 1
    patch = [patch, patch];
end

w = size(mask, 2);
h = size(mask, 1);

pi = patch(1);
pj = patch(2);

% Extract boundary of mask.
dM = conv2(double(mask),[1,1,1;1,-8,1;1,1,1],'same');

dM(dM>0) = 1;
dM(dM<0) = 0;

% Calculate gradients.
% Calculate border confidence.
% Calculate weights.
if ENABLE_LOCAL_WEIGHING
    weight = fspecial('gaussian', G_SIZE, G_SIGMA);
    for i=1:4
        cI = dM.*abs(conv2(II{i},[1,1,1;1,-8,1;1,1,1],'same'));
        if ENABLE_LOCAL_WEIGHING
            wI(i) = { 1./(1+10*conv2(cI, weight, 'same')) };
            %%tt = cI;
            %%tt(isfinite(tt)) = 1;
            %%wI(i) = {tt};
        end
    end
else
    for i=1:4
        cI(i) = { dM.*abs(conv2(II{i},[1,1,1;1,-8,1;1,1,1],'same'))} ;
    end
end

% Allocate I_rec, without setting values (faster)
I_rec(size(mask,1), size(mask,2)) = 0;

% Loop over patches
ni = ceil(h/pi);
nj = ceil(w/pj);

% Debug / stats output.
stats.confidences = [];
stats.avg_nop = 0;

for i=1:ni
    for j=1:nj
        si = (i-1)*pi+1;
        sj = (j-1)*pj+1;
        
        if ENABLE_LOCAL_WEIGHING
            pp = zeros(pi,pj);
            ww = zeros(pi,pj);
            for k = 1:4
                wp = wI{k}(si:si+pi-1, sj:sj+pj-1);
                if all(isfinite(wp))
                    pp = pp + wp.*II{k}(si:si+pi-1, sj:sj+pj-1);
                    ww = ww + wp;
                end
            end
            pp = pp ./ ww;
            I_rec(si:si+pi-1, sj:sj+pj-1) = pp;
            stats.avg_nop = stats.avg_nop + 0;
        else
            conf = zeros(1,4);
            for k=1:4
                conf(k) = sum(sum(cI{k}(si:si+pi-1, sj:sj+pj-1)));
            end

            % Patches with no information have confidence Inf!!!
            [conf, indices] = sort(conf);
            stats.confidences(end+1, :) = conf;

            % Index:
            indWinner = indices(1);

            if conf(1) == 0
                I_rec(si:si+pi-1, sj:sj+pj-1) = II{indWinner}(si:si+pi-1, sj:sj+pj-1);
                stats.avg_nop = stats.avg_nop + 1;
            else
                conf = conf / conf(1);
                indices = indices(conf < CONFIDENCE_THRESHOLD);
                pp = zeros(pi,pj);
                for k = indices
                    pp = pp + II{k}(si:si+pi-1, sj:sj+pj-1);
                end
                pp = pp/length(indices);
                I_rec(si:si+pi-1, sj:sj+pj-1) = pp;
                stats.avg_nop = stats.avg_nop + length(indices);
            end
        end
    end
end

% Code that was used to determine the above threshold
% CC = stats.confidences;
% CC(~isfinite(CC)) = 0;
% BB = CC ./ repmat(CC(:,1),1,4);
% BB(~isfinite(BB)) = 0;
% B = BB(:,2:4); B = sort(B(:));
% figure; plot(B);
% pause;

stats.avg_nop = stats.avg_nop / (ni*nj);
    

end