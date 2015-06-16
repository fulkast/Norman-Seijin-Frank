% Non Greedy sparseCoding

function [Z, U] = epsilonGreedysparseCoding(U, X, M, sigma, rc_min, Im_size)
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

    
    % Initialize the residual with the observation x
    % For the modification with masking make sure that you only take into
    % account the known observations defined by the mask M
    % Initialize z to zero
    % TO BE FILLED
    
    
    m = M(:,1);
    residual = m.*X(:,1);
%     x = residual;
    
%     residual = m.*X;
    x = residual;
    Uloc = U.*(repmat(m,1,size(U,2)));
    rc_max = max((Uloc'*residual).^2);
    
    it = 0;
    residual = M.*X;
while (norm(residual) > sigma*norm(x)) && (rc_max > rc_min)
%         it=it+1;
        % TO BE FILLED 
    it = it +1;
        
    for nn = 1:n  
        m = M(:,nn);
        Uloc = U.*(repmat(m,1,size(U,2)));
        % Select atom with maximum absolute correlation to the residual       
          
%           vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
          
        vec1 = (Uloc'*(residual(:,nn)))./exp(Z(:,nn).^2/(sum(Z(:,nn))+0.00001));
        vec = abs(vec1);
          
        
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        
        % Update the maximum absolute correlation
        % Update coefficient vector z and residual r
        
        residual(:,nn) = residual(:,nn) - vec1(arg)/norm(d)*d;
%         Z(arg,nn)=Z(arg,nn)+exp(-it/1000)*(vec1(arg)/norm(d));
        Z(arg,nn)=Z(arg,nn)+vec1(arg)/norm(d);
        % For the inpainting modification make sure that you only consider
        % the known observations defined by the mask M
        
       % imnew = DictionaryPlot(U*Z(:,nn),16);
%        subplot(1,3,1) 
%         imshow(imresize(reshape(U*Z(:,nn),blocksize,blocksize),16))
%         title(strcat('Current Progress on patch ', num2str(nn)))
%         subplot(1,3,2) 
%         imshow(imresize(reshape(X(:,nn),blocksize,blocksize),16))
%         title('Actual Patch Trying to Fit')
%         subplot(1,3,3) 
%         plot(sort(vec))
%         title('Coherence of Residue to Each Atom (sorted in ascending order)')
%         xlabel('Atoms Sorted')
%         ylabel('Coherence Value')
%         drawnow;
%     %   subplot(1,2,2)
       
      % imshow(DictionaryPlot(U*Z,Im_size,16))
      
      
    end
    if (mod(it,10) == 0)
    
subplot(2,1,2)
%1)
    pcolor(abs(Z))       % Heatmap of Z coefficients showing how "popular" the atoms are 
    title('Absolute Value of Z coefficients')
%2)
%     pcolor((Z))        % Same as above but now shows the sign of the
%     coefficients
%     title('Z coefficients') 
    colorbar
subplot(2,1,1)
%1)
imshow(DictionaryPlot(U*Z,Im_size,blocksize))
%     plot(max(Z,[],2),'b.','MarkerSize',20); xlabel('Atom');ylabel('Maximum Coefficient');title('Performance of Atoms')        % Shows the maximum coefficient of an atom to any and all of the patches 
%2)
%       plot(mean(Z,2),'b.','MarkerSize',20); xlabel('Atom');ylabel('Mean Coefficient');title('Performance of Atoms')
    
%3)    
%     imshow(imresize(reshape(U(:,(max((max(Z,[],2)))==(max(Z,[],2)))),blocksize,blocksize),16));title('Most
%     Popular Atom');
%     % Displays the most popular atom

    
% colormap(parula);
drawnow
% [y,i] = sort(mean(abs(Z),2),'descend');
% [z,j] = sort(sum(abs(residual),1));



    end
    fprintf('Run %d of %d\n',nn,n);
    % Add the calculated coefficient vector z to the overall matrix Z
%     Z(:,nn) = z;


end
