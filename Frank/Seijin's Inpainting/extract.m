function [X,n1,n2]=extract(I,d)
n1=size(I,1);
n2=size(I,2);

Nd1=(floor(size(I,1)/d)+1)*d;
Nd2=(floor(size(I,2)/d)+1)*d;
I1=zeros(Nd1,Nd2);
for  i=1:n1
    for j=1:n2
        I1(i,j)= I(i,j);
    end
end
N1=fix(Nd1/d);
N2=fix(Nd2/d);
X=zeros(d*d,N1*N2);
for i=1:N1
    for j=1:N2
        for k=0:d*d-1
            X(k+1,(i-1)*N2+j)=I1(1+(i-1)*d+fix(k/d),1+(j-1)*d+rem(k,d));
        end
    end
end
end
