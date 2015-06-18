function [Z,newM] = sparseCoding(U, X, M, sigma, rc_min)
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
Y=U'*(X.*M);

for nn = 1:n
    if mod(nn,floor(n/4))==0
        
       	fprintf(' Percentage Complete: %.2f \n',100*nn/n)
        toc
        tic
    end
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    x=X(:,nn);
    z=Z(:,nn);
    m=M(:,nn);
    missing=1; 
    if (missing>0)
    y=Y(:,nn);
    v=U'*(U.*repmat(m,1,l));
  %  sparseness=0;
   rc_max = Inf;
    while (norm(y)>sigma*norm(x)) && (rc_max > rc_min) %&&(sparseness/dim<0.5)
        [foo,k]=max(y.^2,[],1);
  %      if (z(k)==0) sparseness=sparseness+1;
  %      end
        z(k)=z(k)+y(k);
        y=y-y(k)*v(:,k);
         rc_max = abs(y(k));
     end
    end
    Z(:,nn) = z;
end
% Z=Y;
% 
% Xn=repmat(sum(X.^2),l,1);
% id=find(Z.^2<sigma*Xn);
% Z(id)=0;


    
