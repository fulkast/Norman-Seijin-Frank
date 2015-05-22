function [U,Z] = dictionary_learning(X)
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

l = 40;
sigma = 0.1;
iter_num = 15;
init_mode = 'samples';


%% Initialization of Dictionary


if strcmp(init_mode, 'rand')
    
    % Initialize D with random unit length atoms
    
elseif strcmp(init_mode, 'samples')
    
    % Draw uniform samples from data matrix
else
    error('Invalid value for parameter init_mode.')
end



Z = zeros(l,n);
U_new = zeros(size(U));

for i=1:iter_num
    disp(['iteration: ' num2str(i)]);
    
     Z = mp(U, X, sigma, rc_min);   
    
    for a = some_order       
    end
    
    U(:,2:end) = U_new(:,2:end);
    
end
