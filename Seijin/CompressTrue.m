function I_comp = Compress(I)
d=10;
[X1,n1,n2]=extract(I(:,:,1),d);
Nd1=(fix(n1/d)+1)*d;
Nd2=(fix(n2/d)+1)*d;
N1=fix(Nd1/d);
N2=fix(Nd2/d);
C=size(I,3);

%V=zeros(d*d,k,C);
%mu=zeros(d*d,C);
%Z=zeros(k,N1*N2,C);
for c=1:C
    [X1,n1,n2]=extract(I(:,:,c),d);
    [mu1,lambda1,U1]=PCAnalyse(X1);
    k=0;
    sum=0;
    D=diag(lambda1);
    for i=1:d*d
        sum=sum+D(i);
    end
    sum2=0;
    while (sum2/sum<0.99)
       k=k+1;
       sum2=sum2+D(d*d-k+1);
    end
    V1=U1(:,(d*d-k+1):d*d);
    Z1=V1'*(X1-repmat(mu1',1,N1*N2));
    V{c}=V1;
    mu{c}=(mu1);
    Z{c}=Z1;
end
I_comp= {V,mu,Z,n1,n2,C,d,k};



%I_comp.I = I; % this is just a stump to make the evaluation script run, replace it with your code!

end