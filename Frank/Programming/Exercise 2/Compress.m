function I_comp = Compress(I)

imshow(I);

[r, c ,v] = size(I);
d = 1;
k = 1;

X = extract(I,d)
[mu, lambda, U] = PCAanalyse(X);
Z = zeros(k,size(X,2),v);

% Xtil=zeros(d^2,ceil(r/d)*ceil(c/d),v);

for i = 1:v
    
muMat = repmat(mu(:,i),[1 size(X,2)]);

% size(Z(:,:,i))
% size(U(:,end-k+1:end,i))
% size(X(:,:,i))
% size(muMat)
Z(:,:,i) = U(:,end-k+1:end,i)'*(X(:,:,i)-muMat);


% xbartil = U(:,end-k+1:end,i)*Z;


% Xtil(:,:,i) = xbartil+muMat;
end


I_comp.Z = Z;
I_comp.U = U(:,end-k+1:end,:);
I_comp.mu = mu;
% I_comp.compressed = Xtil;
I_comp.aspect = size(I);
I_comp.d = d;

 % this is just a stump to make the evaluation script run, replace it with your code!
