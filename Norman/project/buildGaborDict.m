function U=buildGaborDict(n)
N=90*(n/4+1)^2;
U=zeros(n*n,N);
ct=1;
for mu1=-n/2:4:n/2
    for mu2=-n/2:4:n/2
        for sig=0:2
            for f=0:4
                for tet=0:30:150
                   U(:,ct)=(Gabor(n,mu1,mu2,1/3^sig,f,tet));
                   ct=ct+1;
                end
            end
        end
    end
end
U(:,1)=ones(n*n,1);
U = normc(U);
save('GaborDictionary.mat','U');
end
