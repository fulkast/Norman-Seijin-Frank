function [Kopt,z] = PredictMissingValues(X,Kmin,Kmax,maxIter,T,method) %monte carlo very good
%Now apply k-means
iter=maxIter;
eps=0.001;
varargin=cell(2,1);
varargin{1}='maxiter';
varargin{2}=iter;
%T=2;
Kopt=0;
Kscore=zeros(Kmax-Kmin+1,1);
scoreOpt=realmax;
[dd,N]=size(X);
zopt=zeros(Kmax-Kmin+1,N);
for k=Kmin:Kmax
    score=zeros(T,1);
    optimT=1;
    Zt=zeros(T,N);
    Ztk=zeros(T,k);
    for t=1:T
        [z,Z]=gmm(X, k,varargin{:});
        Zt(t,:)=z(:);
        Ztk(t,:)=sum(Z,2);
    end
  
   
   %Monte Carlo?
   if (method=='MonteCarlo')
   %for t1=1:T-1
   for t=2:T
      error=realmax;
      for tn=1:5*k*log(k)^3
          if (mod(tn,1000)==0)
             % display(tn);
          end
          a=0;
          b=0;
         while(a==b)
              a=fix(k*rand(1))+1;
              b=fix(k*rand(1))+1;
         end
        ida=find(Zt(t,:)==a);
        idb=find(Zt(t,:)==b);
        A=Zt(t,:);
        A(ida)=b;
        A(idb)=a;
        erraux=0;
        for t1=1:t-1
            aux=numel(find(Zt(t1,:)~=A));
            erraux=erraux+aux;
        end
        r=rand(1);
        
        if(r<exp((error-erraux)*(tn/k)^(0.5))&&error~=erraux)
      %  if(error>erraux)
            Zt(t,:)=A;
            error=erraux;
        end

      end
      for t1=1:t-1
          aux=numel(find(Zt(t1,:)~=Zt(t,:)));
          score(t)=score(t)+aux;
          score(t1)=score(t1)+aux;
      end
      Kscore(k-Kmin+1,1)=Kscore(k-Kmin+1,1)+error^(1)/(T*(T-1));
   end
  %end
   [foo,optimT]=min(score,[],1);
   zopt(k-Kmin+1,:)=Zt(optimT,:);
   else
  
     
    %Cuisine method
   A=repmat(Zt,1,1,N)-permute(repmat(Zt,1,1,N),[1,3,2]);
    for t=1:T
       Zt(t,:)=Ztk(t,Zt(t,:)); 
    end
    Zt=repmat(Zt,1,1,N)+permute(repmat(Zt,1,1,N),[1,3,2]);
    Zt=Zt+eps;
    Zt=Zt.^(-1);
   IDX=find(A~=0);
   A(IDX)=1;
   Zt=A.*Zt;
   A=sum(A,1);
   Zt=sum(Zt,1);
   A=Zt.*(T-A);
   Kscore(k-Kmin+1,1)=sum(sum(A))/(T*(T-1));
end   
   
    if (Kscore(k-Kmin+1)<scoreOpt) 
        Kopt=k;
        scoreOpt=Kscore(k-Kmin+1);
    end
end

% update missing values X_pred
Kscore
disp(Kopt)
z(:)=zopt(Kopt-Kmin+1,:);
end
