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
for nn = 1:n
    nn
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    x=X(:,nn);
    z=Z(:,nn);
    m=M(:,nn);
    idx=find(m~=0);
    m(idx)=1;
    nummissing=numel(find(m==0));
    res=x;
    % TO BE FILLED
    ct=0;

    while (sum((U'*(res.*m)).^2)>sigma&&nummissing>0)
        % TO BE FILLED 
        ct=ct+1;
        %sum((res.*m).^2)
    %    y=U'*(res.*m);
        y=U'*(res.*m);       
        % Select atom with maximum absolute correlation to the residual       
        [foo,k]=max(y.^2,[],1);
        % Update the maximum absolute correlation   TODOOOOOO
   
        % y=(U'*res);
         res=res-y(k)*U(:,k);
    
        % Update coefficient vector z and residual z
        z(k)=z(k)+y(k);
        %res=x-U*z; 
        % For the inpainting modification make sure that you only consider
        % the known observations defined by the mask M
        
    end
    ct
    % Add the calculated coefficient vector z to the overall matrix Z
    Z(:,nn) = z;
end
