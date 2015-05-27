function h=haarTrans(N)

% Computes the normalized Haar Wavelet basis
%
% INPUT
% N: The dimension of the basis - it should be a power of 2 to work
% 
% OUTPUT
% H: NxN HAAR transform matrix. H contains the Haar wavelet basis as
%    columns vectors.
%
%
% A Note on 2-dimensional Haar transformations: 
% We can write H'*A for the HAAR transformation of the columns of an 
% input A, and H*A for the inverse transformation of the columns of A 
% (assuming the input image A is N-by-N).
% In this case, i.e. if A is square, the two-dimensional Haar 
% transformation of A can be computed as H'*A*H.

h=zeros(N,N); 
h(1,1:N)=ones(1,N)/sqrt(N);

for k=1:N-1 
  p=fix(log(k)/log(2)); 
  q=k-(2^p); 
  k1=2^p;
  t1=N/k1; 
  k2=2^(p+1);
  t2=N/k2; 
  for i=1:t2 
    h(k+1,i+q*t1) = (2^(p/2))/sqrt(N); 
    h(k+1,i+q*t1+t2) = -(2^(p/2))/sqrt(N); 
  end 
end
h = h';