function [U,Z] = dictionary_learning(X,Isize)
% Implements dictionary learning algorithm, using matching pursuit as the
% sparse coding stage.
%
% INPUTS
% X: (d x n) data matrix (samples as columns)
%
% OUTPUTS
% U: (d x l) dictionary
% Z: sparse coding of X in dictionary U
%
% PARAMETERS
% l: codebook size of dictionary. SHOULD BE A SQUARE NUMBER
% init_mode: initialization of dictionary, either 'rand', 'samples',
% 'recover' to train the dictionary.mat again,'dct' to initialize with
% DCT,'kmean' to apply Kmean on the training set and initialize on the
% centroids
% iter_num: number of update iterations
% sigma: desired maximal residual norm
% rc_min
% project: Set to true if you want to project the training set onto its
% most significant l-dimentional subspace through PCA

%%Fixed parameters
[d ,n] = size(X);
patch=sqrt(d);
ovlp=0;
%% Parameters

l = d/4;%(should be a carree)
Usize=[sqrt(d*l) sqrt(d*l)];
sigma = 0.01;
rc_min = 0.1;
iter_num = 15;
init_mode = 'rand';
project=false;

if project
    mu=mean(X');
    Xx=X-repmat(mu',1,size(X,2));
    Xx=Xx-repmat(mean(Xx),size(X,1),1);
    [U,lambda]=eig(cov(Xx'));
     U=fliplr(U);
     Y=U(:,1:l-1)*U(:,1:l-1)'*(X-repmat((mean(X)),d,1));
     X=Y+repmat((mean(X)),d,1);
end

%% Initialization of Dictionary


if strcmp(init_mode, 'rand')
    U=rand(d,l);
    U=U-0.5;
    U(:,1) = ones(d,1)/d;
    U = normc(U);
    % Initialize D with random unit length atoms
    
elseif strcmp(init_mode, 'samples')
    
    % Draw uniform samples from data matrix
    U = X(:,randsample(1:n,l));
    U=U-repmat((mean(U)),d,1);
    U(:,1) = ones(d,1)/d;
    U = normc(U);
elseif strcmp(init_mode, 'dct')
    U=overDCTdict(d,l);
elseif strcmp(init_mode, 'kmean')
    %if first coordinate is negative, inverse the vector, to avoid
    %redundant centroids in terms of linear combination
    Y=normc(X-repmat((mean(X)),d,1));
    id=find(Y(1,:)<0);
    Y(:,id)=-Y(:,id);
    %Kopt applies kmean 10 times and chooses the best clustering
     [foo,z] = Kopt(Y,l-1,l-1,25,10,'MonteCarlo');
     Zz = zeros(l-1, n);
     for i=1:n
         Zz(z(1,i),i)=1;
     end
    U = Y*Zz';
    U = U./repmat(sum(Zz,2)'+0.01,d,1);
    V = ones(d,l);
    V(:,2:l)=U(:,:);
    U = normc(V);

elseif strcmp(init_mode, 'recover')
    t=load('dictionary.mat');
    U=t.U;
else
    error('Invalid value for parameter init_mode.')
end
subplot(1,2,1)
imshow(my_col2im(X,patch,Isize,ovlp,[0 0]))  
subplot(1,2,2)
imshow(0.5+my_col2im(0.5*(U-repmat((mean(U)),d,1))./repmat(sqrt(var(U)),d,1),patch,Usize,0,[0 0]));
figure
Z = zeros(l,n);

U_new = zeros(size(U)); 
for i=1:iter_num

    disp(['iteration: ' num2str(i)]);
    
    Z = sparseCoding(U, X, ones(size(X)) , sigma, rc_min);   
    subplot(1,2,1)
    imshow(my_col2im(U*Z,16,Isize,ovlp,[0 0]))  
    subplot(1,2,2)
imshow(0.5+my_col2im(0.5*(U-repmat((mean(U)),d,1))./repmat(sqrt(var(U)),d,1),16,Usize,0,[0 0]));
    for a = 2:l 
    
    Nset = find(Z(a,:)~=0);

    if ~isempty(Nset)
    U(:,a) = zeros(size(U,1),1);
    R = X(:,Nset) - U*Z(:,Nset);
    
    g = Z(a,Nset)';
    for it = 3
    h = R*g/norm(R*g);
    g = R'*h;
    U_new(:,a) = h;
    Z(a,Nset) = g';
    end
    U(:,a) = U_new(:,a);  
    end 
    
    end
    
figure
  
title(num2str(norm(X-U*Z)))
save('dictionary.mat','U');
end
imshow(my_col2im(U*Z,patch,Isize,ovlp,[0 0]))      

save('dictionary.mat','U');
end
