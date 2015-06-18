% This one goes incrementally from higher confidence to lower

function I_out = EM4All( X, M)

masking_quality_cutoff = .5;

sigma = 0.09; %stopping criterion for sparsecoding
rc_min = 0.0003; %stopping criterion for sparsecoding 


d = 16;

    [X, Im_size]= extract(X,d);
    
    [M, ~]= extract(M,d);
       
    Yobs = X;
    Yobs(M==0) = 0;
    % Get dictionary 

    U = buildDictionary(d^2,0);  %    0 second argument for CUSTOM Dictionary
    
blocksize = size(X,1)^.5;
l = size(U,2);
n = size(X,2);
Z = zeros(l,n);


% Loop over all observations in the columns of X

[maskcount, order] = sort(sum(M,1),'descend');


seq = order;
seq(maskcount/blocksize/blocksize >= masking_quality_cutoff) = [];

% order(maskcount/blocksize/blocksize == 1) = [];

for nn = order
    % TO BE FILLED
    
    Yobsloc = Yobs(:,nn);
    converged = 0;
    
    %% Initial sparse coding used to feed EM
    m = M(:,nn);
    residual = m.*X(:,nn);
    Uloc = U.*(repmat(m,1,size(U,2)));
    
    Z = SCoding(Uloc,Z,residual,sigma,rc_min,nn);
    
    %% continue with EM
    m = M(:,nn);
    residual = Yobsloc +(1-m).*(U*Z(:,nn));
    Z2 = EM(m,residual,Yobsloc,U,converged,sigma,rc_min,l);    
    Z(:,nn) = Z2;
    
    %% work on neighbours
    

end

Yobs = Yobs + (1-M).*(U*Z);

I_out = DictionaryPlot(Yobs,Im_size,blocksize);

end


