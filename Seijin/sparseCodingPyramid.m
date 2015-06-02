function Z = sparseCodingPyramid(U, X, M, sigma, rc_min)


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
    if mod(nn,1000)==0
        nn
    end

    x=X(:,nn);
    z=Z(:,nn);

      m=M(:,nn);
      idx=find(m<=0);
      m(idx)=0;

    % TO BE FILLED
    ct=0;
    elemt=elemt+1;
    y=Y(:,nn);
    v=U'*(U.*repmat(m,1,size(U,2)));
    sparseness=0;
    while (sum(y.^2)>sigma*norm(x))%&&(sparseness/dim<0.5)
        % TO BE FILLED 
        ct=ct+1;   
        [foo,k]=max(y.^2,[],1);
        if (z(k)==0) sparseness=sparseness+1;
        end
        z(k)=z(k)+y(k);
        y=y-y(k)*v(:,k);
        
    end
 
    % Add the calculated coefficient vector z to the overall matrix Z
    Z(:,nn) = z;
end

