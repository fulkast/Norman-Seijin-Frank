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

l = size(U,2);
n = size(X,2);
Z = zeros(l,n);
% Loop over all observations in the columns of X
tic
for nn = 1:n
    if mod(nn,floor(n/4))==0
        
       	fprintf(' Percentage Complete: %.2f \n',100*nn/n)
        toc
        tic
    end
    % Initialize the residual with the observation x. Only take into  
    % account the known observations that are not masked out by m
    m = M(:,nn);
    miss=numel(find(m==0));
    idx=find(m<=0);
    m(idx)=0;
    x = X(m~=0,nn); % all entries of column nn, except the masked ones.
    u = U(m~=0,:);  % all rows of U, except except the masked ones.
    r = x;
    z = zeros(size(U,2),1);
    rc_max = Inf;

    while (norm(r) > sigma*norm(x)) && (rc_max > rc_min)&&(miss>0)
        
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
end
toc
