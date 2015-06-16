
%step=5;
%ini=0;
Y=zeros(3,1);
for i=1:10
    r=5*rand(1);
    nu=10+fix(30*rand(1));
    X=r*rand(3,nu);
    %ini=ini+step;
    ini=30*rand(3,1);
    elong=2*rand(3,1);
    X1=repmat(elong,1,nu).*X+repmat(ini,1,nu);
    Y=[Y X1];
end
ct=0;
N=1;
K=zeros(N,1);
z=zeros(size(Y,2),1);
for i=1:N
[K(i),z]=Kopt(Y,6,13,100,5,'MonteCarlo');
end
mean1=mean(K)
var1=var(K)
k=round(mean1)
%[K(i),z]=Kopt(Y,k,k,1000,30,'MonteCarlo');
varargin={'maxiter',1000};
Y=Y';
scatter3(Y(:,1),Y(:,2),Y(:,3),10,z,'filled');