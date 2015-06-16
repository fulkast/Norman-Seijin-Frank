function [I_rec, stats] = patch_selector(I1, I2, I3, I4, mask, patch)
% stats optionally informs caller about some statistics
%   stats.avg_nop:      average number of patches used for blending
%   stats.confidences:  confidence data per patch

CONFIDENCE_THRESHOLD = 1.5;

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
cI1 = abs(conv2(double(I1),[1,1,1;1,-8,1;1,1,1],'same'));
cI2 = abs(conv2(double(I2),[1,1,1;1,-8,1;1,1,1],'same'));
cI3 = abs(conv2(double(I3),[1,1,1;1,-8,1;1,1,1],'same'));
cI4 = abs(conv2(double(I4),[1,1,1;1,-8,1;1,1,1],'same'));

% Calculate border confidence
cI1 = dM.*cI1;
cI2 = dM.*cI2;
cI3 = dM.*cI3;
cI4 = dM.*cI4;

maxConf = [max(cI1(:)), max(cI2(:)), max(cI3(:)), max(cI4(:))];

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
        conf = zeros(1,4);
        conf(1) = sum(sum(cI1(si:si+pi-1, sj:sj+pj-1)));
        conf(2) = sum(sum(cI2(si:si+pi-1, sj:sj+pj-1)));
        conf(3) = sum(sum(cI3(si:si+pi-1, sj:sj+pj-1)));
        conf(4) = sum(sum(cI4(si:si+pi-1, sj:sj+pj-1)));
        
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