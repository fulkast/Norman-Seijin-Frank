function [Z,newM] = sparseCoding(U, X, M,MS, sigma, rc_min)
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
dim=size(X,1);
Z = zeros(l,n);
% Loop over all observations in the columns of X
newM=ones(size(M));

elemt=0;
INFO=zeros(n);
Y=U'*(X.*M);

for nn = 1:n
    if mod(nn,500)==0
        nn
    end
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    x=X(:,nn);
    z=Z(:,nn);
   ms=MS(:,nn);    
   info=numel(find(ms~=0))/dim;
    m=M(:,nn);
    missing=1;
   missing=dim-numel(find(m~=0));
      idx=find(m<=0);
      m(idx)=0;
   infoMbis=numel(find(m~=0));
   moy=sum(x.*m)/infoMbis;
   vari=sum((x-moy).^2.*m)/infoMbis;

    

    info=(sqrt(vari))/(info)^2;
    INFO(nn)=info;
    % TO BE FILLED
    ct=0;
    if (missing>0)%&&numel(find(ms~=0))/dim>0.75)
    elemt=elemt+1;
    y=Y(:,nn);
    v=U'*(U.*repmat(m,1,size(U,2)));
    while (sum(y.^2)>sigma)
        % TO BE FILLED 
        ct=ct+1;
        %sum((res.*m).^2)
        %y=U'*(res.*m);       
        % Select atom with maximum absolute correlation to the residual       
        [foo,k]=max(y.^2,[],1);
        % Update the maximum absolute correlation   TODOOOOOO
         %res=res-y(k)*U(:,k);
        % Update coefficient vector z and residual z
        z(k)=z(k)+y(k);
        y=y-y(k)*v(:,k);
        %res=x-U*z; 
        % For the inpainting modification make sure that you only consider
        % the known observations defined by the mask M
        
     end
    else
        newM(:,nn)=m;
    end
 
    % Add the calculated coefficient vector z to the overall matrix Z
    Z(:,nn) = z;
end
% 
for i=1:floor(elemt*0.003*0)
    [foo,k]=max(INFO,[],1);
    newM(:,k)=M(:,k);
    INFO(k)=0;
end
    
