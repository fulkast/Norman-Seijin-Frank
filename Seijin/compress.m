function Ic=compress(I,n)
Ic=I;
if (n~=1)
N1=floor(size(I,1)/n);
N2=floor(size(I,2)/n);
Ic=zeros(N1,N2,size(I,3));
for i=1:N1
    for j=1:N2
        Ic(i,j,:)=mean(mean(I((i-1)*n+1:n*i,(j-1)*n+1:n*j,:)));
    end
end
end