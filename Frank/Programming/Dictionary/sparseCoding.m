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
    
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    m = M(:,nn);
    residual = m.*X(:,nn);
    x = residual;
    
    Uloc = U.*(repmat(m,1,size(U,2)));
    
    rc_max = max((Uloc'*residual).^2);
    % TO BE FILLED
    
    while (norm(residual) > sigma*norm(x)) && (rc_max > rc_min)
        
        % TO BE FILLED 
        
        % Select atom with maximum absolute correlation to the residual       
          
          vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
          vec = abs(vec1);
          
        
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        
        % Update the maximum absolute correlation
        % Update coefficient vector z and residual r
        
        residual = residual - vec1(arg)/norm(d)*d;
        Z(arg,nn)=Z(arg,nn)+vec1(arg)/norm(d);
        % For the inpainting modification make sure that you only consider
        % the known observations defined by the mask M
%         
%         for l = 1: length(Uloc)
%             R = x - Uloc*Z(:,nn);
%             g = Z(l,nn);
%             h = R*g/norm(R*g);
%             g = R'*h;
%             Uloc(:,l) = h;
%             Z(l,nn) = g';
%         end
        
        
    end
    
    
    fprintf('Run %d of %d\n',nn,n);
    % Add the calculated coefficient vector z to the overall matrix Z
%     Z(:,nn) = z;
end
