function Ibig=expand(I,n)

N1=n*size(I,1);
N2=n*size(I,2);
Ibig=zeros(N1,N2,size(I,3));
for i=1:N1
    for j=1:N2
        Ibig(i,j,:)=I(floor((i-1)/n+1),floor((j-1)/n+1),:);
    end
end
        
