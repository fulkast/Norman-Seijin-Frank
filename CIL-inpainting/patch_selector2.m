function [I_rec, stats] = patch_selector2(I1, I2, I3, I4, mask, patch)
% stats optionally informs caller about some statistics
%   stats.avg_nop:      average number of patches used for blending
%   stats.confidences:  confidence data per patch
%
% WATCH OUT: I1, I2, I3, I4 should contain NaNs for those patches that are
% not complete (sides of the images). This makes the interface dirty, but
% since this is just a test...

CONFIDENCE_THRESHOLD = 1.5;

II = {I1, I2, I3, I4};

if numel(patch) == 1
    patch = [patch, patch];
end

w = size(mask, 2);
h = size(mask, 1);

pi = patch(1);
pj = patch(2);

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
        
        
        
        % Extract (inner) boundary of mask.
        M = mask(si:si+pi-1, sj:sj+pj-1);
        dM = conv2(double(M),[1,1,1;1,-8,1;1,1,1],'same');
        dM(dM>0) = 1;
        dM(dM<0) = 0;
        dMI = (dM~=0);
        
        [gMx, gMy] = gradient(double(M));
        gM = [gMx(dMI), gMy(dMI)];
        gM = normr(gM);
        
        % Calculate border confidence
        C = [];
        isInf = false(1,4);
        for k=1:4
            p = II{k}(si:si+pi-1, sj:sj+pj-1);
            isInf(k) = any(~isfinite(p(:)));
            [gIx, gIy] = gradient(p);
            gI = [gIx(dMI), gIy(dMI)];
            c = sum(gM.*gI,2);
            C(:,k) = abs(c);
        end
        
        conf = sum(C,1);
        conf(isInf) = Inf;
        
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
                p = II{k}(si:si+pi-1, sj:sj+pj-1);
                if (any(~isfinite(p)))
                    disp here;
                end
                pp = pp + II{k}(si:si+pi-1, sj:sj+pj-1);
            end
            pp = pp/length(indices);
            if (any(~isfinite(pp)))
                disp here;
            end
            I_rec(si:si+pi-1, sj:sj+pj-1) = pp;
            stats.avg_nop = stats.avg_nop + length(indices);
        end
        pp = I_rec(si:si+pi-1, sj:sj+pj-1);
        if (any(~isfinite(pp)))
            disp here;
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