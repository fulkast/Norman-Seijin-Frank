function I_rec = Decompress(I_comp)
V=I_comp{1};
mu=I_comp{2};
Z=I_comp{3};
n1=I_comp{4};
n2=I_comp{5};
C=I_comp{6};
d=I_comp{7};
k=I_comp{8};

% Your decompression code goes here!
Nd1=(fix(n1/d)+1)*d;
Nd2=(fix(n2/d)+1)*d;
N1=fix(Nd1/d);
N2=fix(Nd2/d);
Y=zeros(d*d,N1*N2,C);   
for c=1:C
    Y(:,:,c)=V{c}*Z{c}+repmat(mu{c}',1,N1*N2);
end
    Ix=(zeros(N1,N2,C));
for i=1:N1
    for j=1:N2
        for l=0:d*d-1
            for c=1:C
              Ix(1+(i-1)*d+fix(l/d),1+(j-1)*d+rem(l,d),c)=Y(l+1,(i-1)*N2+j,c);    
            end
        end
    end
end
I_rec =(Ix(1:n1,1:n2,:));%im2uint8(Ix(1:n1,1:n2,:)); % this is just a stump to make the evaluation script run, replace it with your code!