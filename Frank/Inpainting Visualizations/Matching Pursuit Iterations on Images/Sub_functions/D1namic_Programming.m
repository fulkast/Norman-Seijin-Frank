% Dynamic Programming for Sparse Coding

function Z = D1namic_Programming(U, X, M, sigma, rc_min, Im_size)
% 
% Perform sparse coding using a modified matching pursuit tailored to the 
% inpainting problem with residual stopping criterion.
%
% INPUT
% U: (d x l) unit norm atoms
% X: (d x n) observations
% M: (d x n) mask denoting which observations are unknown
% sigma: residual error stopping criterion, normalized by signal norm
% rc_min: minimal residual correlation before stopping
%
% OUTPUT
% Z: MP coding
%


l = size(U,2);
n = size(X,2);
blocksize = size(X,1)^.5;

Z = zeros(l,n);
% Loop over all observations in the columns of X
for nn = 1:n
    
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    m = M(:,nn);
    residual = m.*X(:,nn);
    x = residual;
    
    Uloc = U.*(repmat(m,1,size(U,2)));
    
    rc_max = max((Uloc'*residual).^2);
    % TO BE FILLED
    it = 0;
    conv = 0;
    while (norm(residual) > sigma*norm(x)) && (rc_max > rc_min)
        
        [arg, innerproduct,secondarg,secondinner,conv,rc_max] = S4condStepisGreedy(Uloc,residual,conv);
        d = Uloc(:,arg(1));
        residual = residual - innerproduct/norm((d))*d;
        Z(arg(1),nn)=Z(arg(1),nn)+innerproduct/norm(d);
        
        
        if ~conv
        d = Uloc(:,secondarg(1));
        residual = residual - secondinner/norm((d))*d;
        Z(secondarg(1),nn)=Z(secondarg(1),nn)+secondinner/norm(d);
        end

        it = it +1;
        
%         ProgPlotter(nn,U,Z,X,vec,blocksize)
%         fprintf('Iteration number %.2d\n',it)
%         imshow(imresize(reshape(U*Z(:,nn),blocksize,blocksize),16))
        
% 
%         subplot(1,4,1) 
% imshow(imresize(reshape...
% (U(:,arg)./max(U(:,arg)),blocksize,blocksize),16))
% % Note the 16 above is just for magnifying the image and is not there
% % because of the specific patch size
%         title('Step 1')
%         subplot(1,4,2) 
% imshow(imresize(reshape...
% (U(:,secondarg)./mean(max(U(:,secondarg))),blocksize,blocksize),16))
% % Note the 16 above is just for magnifying the image and is not there
% % because of the specific patch size
%         title('Step 2')
%         subplot(1,4,3) 
% imshow(imresize(reshape(X(:,nn),blocksize,blocksize),16))
%         title('Actual Patch Trying to Fit')
%         subplot(1,4,4) 
% imshow(imresize(reshape...
% (U*Z(:,nn),blocksize,blocksize),16))
%         title('Current Progress')
% %         drawnow;


%         norm(residual)

    end
    
    nn
    fprintf('Took %d iterations\n',it);
    % Add the calculated coefficient vector z to the overall matrix Z
%     Z(:,nn) = z;
end
