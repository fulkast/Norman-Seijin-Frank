function [U,Z,U_new] = dictionary_learning(X)
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

l = 10;
sigma = .01;
rc_min = .01;
% iter_num = 20;
iter_num = 15;
init_mode = 'samples';


%% Initialization of Dictionary


if strcmp(init_mode, 'rand')
    
    % Initialize D with random unit length atoms
    
elseif strcmp(init_mode, 'samples')
    
    % Draw uniform samples from data matrix
    U = X(:,randsample(1:n,l));
    U(:,1) = ones(d,1)/d;
    U = U./(repmat(sqrt(sum(U.^2,1)),d,1));
    
else
    error('Invalid value for parameter init_mode.')
end

% U = buildDictionary(16*16);

Z = zeros(l,n);
U_new = zeros(size(U));
 
for i=1:iter_num
    disp(['iteration: ' num2str(i)]);
    
    Z = sparseCoding(U, X, ones(size(X)) , sigma, rc_min);   
    
    for a = 1:l 
    
    Nset = find(Z(a,:)~=0);
    U(:,l) = zeros(size(U,1),1);
    if ~isempty(Nset)
    R = X(:,Nset) - U*Z(:,Nset);
    
    g = Z(a,Nset)';
    for it = 1
    h = R*g/norm(R*g);
    g = R'*h;
    U_new(:,a) = h;
    Z(a,Nset) = g';
    end
    end
    U(:,1:end) = U_new(:,1:end);    
    
    end
    

figure
imshow(my_col2im(U*Z,16,[512 512]))    
title(num2str(norm(X-U*Z)))
end
    
% U = U_new;
% U(:,2:end) = U_new(:,2:end);
    
end
