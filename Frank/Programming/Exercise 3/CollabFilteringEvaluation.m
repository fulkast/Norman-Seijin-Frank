% Evaluation script for the collaborative filtering problem. This is
% essentially the same script as the evaluation server will run for your 
% submission.
%
% It first loads an input data matrix from Data.mat, and splits the known 
% values (ratings from 1 to 5 stars) into training and testing sets. 
% It then passes the training data to your custom implementation of the 
% function 
%      PredictMissingValues.m
% in order to train, and obtain predictions for all unknown entries.
% Finally, the script compares these prediced entries to the test set,
% and computes the mean squared error (MSE).

% Setup
rand('seed', 1);  % fix random seed for reproducibility

% Constants
filename = 'Data.mat';
prc_trn = 0.5;  % percentage of training data
nil = 0;  % missing value indicator. could also use NaN or any other value instead

% Load data
L = load(filename);
X = L.X;

idx = find(X ~= nil); % indices of existing (non-zero) ratings
n = numel(idx);

% Split into training and testing index sets
n_trn = round(n*prc_trn);
rp = randperm(n);
idx_trn = idx(rp(1:n_trn));
idx_tst = idx(rp(n_trn+1:end));

% Build training and testing matrices
X_trn = ones(size(X))*nil;
X_trn(idx_trn) = X(idx_trn);  % add known training values

X_tst = ones(size(X))*nil;
X_tst(idx_tst) = X(idx_tst);  % add known testing values


% Predict the missing values here!
X_pred = PredictMissingValues(X_trn, nil);


% Compute MSE on the test set
rmse = sqrt(mean((X_tst(X_tst ~= nil) - X_pred(X_tst ~= nil)).^2));  % error on known test values

disp(['Root of Mean-squared error on test set: ' num2str(rmse)]);
