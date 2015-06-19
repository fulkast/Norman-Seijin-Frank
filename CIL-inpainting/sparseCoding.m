function Z = sparseCoding(U, X, M, sigma, rc_min)
% 
% Perform sparse coding using a modified matching pursuit tailored to the 
% inpainting problem with residual stopping criterion.
%
% INPUT
% U: (d x l) unit norm atoms
% X: (d x n) observations
% M: (d x n) mask denoting which observations are unknown
% sigma: residual error stopping criterion, normalized by signal norm
% rc_min: minimal residual correlation before stopping
%
% OUTPUT
% Z: MP coding
%
PROPAGATION = false;

l = size(U,2);
n = size(X,2);
Z = zeros(l,n);
% Loop over all observations in the columns of X

% Set up patch propagation, if enabled.
if PROPAGATION
    blocksize = size(X,1)^.5;
    [maskcount, order] = sort(sum(M,1),'descend');
    masking_quality_cutoff = .5;

    seq = order;
    seq(maskcount/blocksize/blocksize >= masking_quality_cutoff) = [];
else
    order = 1:n;
end

for nn = order
        
    % Initialize the residual with the observation x. Only take into  
    % account the known observations that are not masked out by m
    m = M(:,nn);
    x = X(m~=0,nn); % all entries of column nn, except the masked ones.
    u = U(m~=0,:);  % all rows of U, except except the masked ones.
    r = x;
    z = zeros(size(U,2),1);
    rc_max = Inf;
    
    while (norm(r) > sigma*norm(x)) && (rc_max > rc_min)
        
        % Select atom with maximum absolute correlation to the residual
        % - We get the correlation (the scalar product) by U'*x = (x'*U)'
        % - find max entry.
        corr = r'*u;
        [~,zi] = max(abs(corr));
        
        % Update coefficient vector z and residual r
        z(zi) = z(zi)+corr(zi);
        r = r - corr(zi)*u(:,zi);
        
        % Update the maximum absolute correlation
        rc_max = abs(corr(zi));
        
        % Only consider the known observations defined by the mask M
    end
    % Add the calculated coefficient vector z to the overall matrix Z
    Z(:,nn) = z;
    
    if PROPAGATION
        sigma__ = 0.005; %stopping criterion for sparsecoding
        rc_min__ = 0.3; %stopping criterion for sparsecoding
        X(:,nn) = U*Z(:,nn);
        sz = [512, 512];
        [X,M,seq] = stencil(X,M,U,nn,blocksize,sz,seq,masking_quality_cutoff,sigma__,rc_min__);
        seq(seq==nn) = [];
        figure(1); imshow(my_col2im(X, sqrt(size(X,1)), sz, 8, [0,0], true));% pause(0.1)
        %     subplot(1,2,1)
        %     imshow(DictionaryPlot(X,[32 32],16))
        %     drawnow;
        %     subplot(1,2,2)
        %    imshow(imresize(reshape(X(:,nn),16,16),32));pause
        %     fprintf('Took %d iterations\n',it);
        %     seq = order;
    end
end
