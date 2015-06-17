% This one goes incrementally from higher confidence to lower

function Z = C0nfidence_first(U, X, M, sigma, rc_min, Im_size)
blocksize = size(X,1)^.5;

U1=haarTrans(blocksize*blocksize);


l = size(U,2);
n = size(X,2);
Z = zeros(l,n);


% Loop over all observations in the columns of X

[maskcount, order] = sort(sum(M,1),'descend');

masking_quality_cutoff = 0;

seq = order;
seq(maskcount/blocksize/blocksize >= masking_quality_cutoff) = [];

% order(maskcount/blocksize/blocksize == 1) = [];

for nn = order
% for nn = 1:n
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
    
%     X(:,nn) = U*Z(:,nn);
    
%     imshow(DictionaryPlot(X,[32 32],16)) ; ...pause;
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



end


