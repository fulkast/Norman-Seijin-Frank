function Ans=filterImage(M,condition) %condition=3*3 0,1 matrix, tells the min. cdt so that Mij is conserved at 0. (Mij = 0 if all coord of cdt which is 1 is 0 )
n1=size(M,1);
n2=size(M,2);
m1=size(condition,1);
m2=size(condition,2);
Ans=zeros(n1,n2);
aux=zeros(n1+m1-1,n2+m2-1);
M=im2double(M);
for i=1:m1
    for j=1:m2
        aux(i:(n1+i-1),j:(n2+j-1))=aux(i:(n1+i-1),j:(n2+j-1))+condition(m1+1-i,m2+1-j).*M(:,:);
    end
end

Ans(:,:)=aux(2:n1+1,2:n2+1);





end