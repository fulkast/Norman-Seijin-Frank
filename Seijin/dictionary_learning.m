function [U,Z,U_new] = dictionary_learning(X,Isize)
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
% l: codebook size of dictionary
% init_mode: initialization of dictionary, either 'rand' or 'samples'
% iter_num: number of update iterations
% sigma: desired maximal residual norm

%% Parameters

[d ,n] = size(X);
l = 2*d;
Usize=[d 2*d];
sigma = 0.3;
rate=0.7;
rc_min = 0.01;
% iter_num = 20;
iter_num = 15;
init_mode = 'rand';
ovlp=0;

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
    U(:,1) = ones(d,1)/d;
    
    U = normc(U);
elseif strcmp(init_mode, 'dct')
    U=overDCTdict(d,l);
else
    error('Invalid value for parameter init_mode.')
end

% U = buildDictionary(16*16);

Z = zeros(l,n);

U_new = zeros(size(U)); 
for i=1:iter_num

    disp(['iteration: ' num2str(i)]);
    
    Z = sparseCoding(U, X, ones(size(X)) , sigma, rc_min);   
    subplot(1,2,1)
    imshow(my_col2im(U*Z,16,Isize,ovlp,[0 0]))  
    subplot(1,2,2)
    imshow(0.5+my_col2im(U./repmat(sqrt(std(U)),d,1),16,Usize,0,[0 0]))
    sigma=max(sigma*rate,0.01);
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
end
imshow(my_col2im(U*Z,16,Isize,ovlp,[0 0]))      
% U = U_new;
% U(:,2:end) = U_new(:,2:end);
save('dictionary.mat','U');
end