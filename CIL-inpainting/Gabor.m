function M=Gabor(n,mu1,mu2,sig,f,tet)
M=zeros(n,n);
f=f/pi;
sig=sig*n^2/2;
tet=tet/180*pi;
for i=1:n
    for j=1:n
        i=i-n/2;
        j=j-n/2;
        M(i+n/2,j+n/2)=exp((-(i-mu1)^2-(j-mu2)^2)/sig)*cos(f*(i*cos(tet)+j*sin(tet)));
        i=i+n/2;
        j=j+n/2;
    end
end
M=reshape(M,[n*n,1]);
