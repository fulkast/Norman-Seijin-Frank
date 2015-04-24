function X_pred = PredictMissingValues(X, nil)
% Predict all missing entries (ratings) in matrix X based on the known entries. 
%
% In other words, the known entries in X are your training given, and
% are passed to the function here as the input.
%
% Missing values in X are denoted by the special constant value nil.
% 
% Competition Setup: This is the function you submit to the competition
% system. The evaluation server will use almost the same script as 
%         CollabFilteringEvaluation.m
% but using a different train/test dataset.

[num_users, num_items] = size(X);
idx = find(X ~= nil);

% your collaborative filtering code here!
Xnew = zeros(size(X));
Xnew(idx)=1;
sam = sum(Xnew,2);

% X_pred = repmat(mean(X(X~=0),2),1,size(X,2));
X_pred = repmat(sum(X,2)./sam,1,size(X,2));
X_pred(sam==0,:)=0;
X_pred(idx) = X(idx);

[U,D,V] = svd(X_pred);
Udash = (U(:,1:size(D,2)) .* repmat((diag(D).^0.5)',size(U,1),1));
Vdash = repmat((diag(D).^0.5)',size(V,1),1).*V;

% error = mean2(Udash(:,2:end)*Vdash(:,2:end)'-X_pred);
% error = mean2(Udash*Vdash'-X_pred);
X_pred = Udash(:,[1 5 20])*Vdash(:,[1 5 20])';
% display(error)

% X_pred = X;
% X_pred(X_pred == nil) = 4; % dummy prediction, 4 stars always