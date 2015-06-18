% This one goes incrementally from higher confidence to lower

function I_out = C0nfidence_first( X, M)

sigma = 0.005; %stopping criterion for sparsecoding
rc_min = 0.3; %stopping criterion for sparsecoding 


d = 16;

    [X, Im_size]= extract(X,d);
    [M, ~]= extract(M,d);
    
    % Get dictionary 
    
    I_out = X;
    oldMask = M;
    U = buildDictionary(d^2,0);  %    0 second argument for CUSTOM Dictionary
    
    blocksize = size(X,1)^.5;
% U1=haarTrans(blocksize*blocksize);


l = size(U,2);
n = size(X,2);
Z = zeros(l,n);


% Loop over all observations in the columns of X

[maskcount, order] = sort(sum(M,1),'descend');
masking_quality_cutoff = .5;

seq = order;
seq(maskcount/blocksize/blocksize >= masking_quality_cutoff) = [];

% order(maskcount/blocksize/blocksize == 1) = [];

for nn = order
% for nn = 1:n

    Yobs = X;
    Yobs(M==0) = 0;
    m = M(:,nn);
    residual = m.*X(:,nn);
    x = residual;
    
    Uloc = U.*(repmat(m,1,size(U,2)));
    
    rc_max = max((Uloc'*residual).^2);
    % TO BE FILLED
    it = 0;
% if nn == 1023 NeighbourImprovement_check(M(:,nn),X(:,nn),16);pause; end
    while (norm(residual) > sigma*norm(x)) && (rc_max > rc_min)


        vec1 = (Uloc'*(residual));
        vec = abs(vec1);
        
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Z(arg,nn)=Z(arg,nn)+vec1(arg)/norm(d);
        it = it +1;

    end
    
    %% E-M SECTION
    converged = 0;
    Yobsloc = Yobs(:,nn);
    m = M(:,nn);
    while ~converged
    
    
    residual = Yobsloc +(1-m).*(U*Z(:,nn));
    
%     residual = m.*X(:,nn);
    x = residual;
%     Uloc = U.*(repmat(m,1,size(U,2)));
    Uloc = U;
    rc_max = max((Uloc'*residual).^2);
    
%    imshow(imresize(reshape(residual,16,16),16));pause
   Z2 = zeros(l,1);
    % Maximization step
    while (norm(residual) > 0.1*sigma*norm(x)) && (rc_max > rc_min)
        vec1 = (Uloc'*(residual));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Z2(arg)=Z2(arg)+vec1(arg)/norm(d);
    
    end
        
        
        it = it +1;
        if norm((Yobsloc+(1-m).*(U*Z2))-x)/norm(x) < sigma
            converged = 1;
        end
           
        Z(:,nn) = Z2;
    end
    
%     X(:,nn) = U*Z(:,nn);
    
    imshow(DictionaryPlot(X,[32 32],16)) ; ...pause;
    M(:,nn) = ones(blocksize*blocksize,1);
    X(:,nn) = U*Z(:,nn);
    
    [X,M,seq] = stencil(X,M,U,nn,blocksize,Im_size,seq,masking_quality_cutoff,sigma,rc_min);

    
    seq(seq==nn) = [];
%     subplot(1,2,1)
%     imshow(DictionaryPlot(X,[32 32],16))
%     drawnow;
%     subplot(1,2,2)
%    imshow(imresize(reshape(X(:,nn),16,16),32));pause
%     fprintf('Took %d iterations\n',it);
%     seq = order;
end

I_out = I_out + (1-oldMask).*(U*Z);
I_out = DictionaryPlot(I_out,Im_size,blocksize);

end


