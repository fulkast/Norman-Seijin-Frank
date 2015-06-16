%stencil

function [X,M,seq] =  stencil(X,M,U,nn,blocksize,Im_size,seq,masking_quality_cutoff,neibounorm,rc_min)
Mask_retention = 0; % 0 for all ones
intrusion = floor(1);
last_el = blocksize^2;
blockcutoff = masking_quality_cutoff * intrusion * blocksize;
rc_max = 1;

    topsection = 1:blocksize;
    topsection = topsection-1;
    topsection = topsection * blocksize;
    topsection = repmat(topsection,(blocksize-intrusion),1);
    topsection = topsection + repmat([1:blocksize-intrusion]',1,blocksize);
    topbase = reshape(topsection,numel(topsection),1);
    topsection = 1:blocksize;
    topsection = topsection-1;
    topsection = topsection * blocksize;
    topsection = repmat(topsection,(intrusion),1);
    topsection = topsection + repmat([1:intrusion]',1,blocksize);
    topsection = reshape(topsection,numel(topsection),1);
    topintrusion = topsection;
    bottomintrusion = reshape(blocksize*blocksize - (repmat(fliplr(1:blocksize)-1,intrusion,1)*blocksize + repmat(fliplr(0:intrusion-1)',1,blocksize)),blocksize*intrusion,1);
    bottombase = blocksize*blocksize - reshape((repmat(fliplr(1:blocksize)-1,blocksize-intrusion,1)*blocksize + repmat(fliplr(0:blocksize-intrusion-1)',1,blocksize)),blocksize*(blocksize-intrusion),1);


up = nn-1;
if ~sum(up==seq)
    up = 0;
elseif sum(M(bottomintrusion,up)) > blockcutoff || sum(M(topintrusion,nn)) < sum(M(bottomintrusion,up))
    up = 0;
end

down = nn+1;
if ~sum(down==seq)
    down = 0;
elseif sum(M(topintrusion,down)) > blockcutoff || sum(M(bottomintrusion,nn)) < sum(M(topintrusion,down))
    down = 0;
    
end
left = nn - Im_size(1); 
if ~sum(left==seq) 
    left = 0; 
elseif sum(M(end-intrusion*blocksize+1:end,left)) > blockcutoff ||sum(M(1:intrusion*blocksize,nn)) < sum(M(end-intrusion*blocksize+1:end,left))
    left = 0;
end
right = nn + Im_size(1);
if ~sum(right==seq) 
    right = 0;
elseif  sum(M(1:intrusion*blocksize,right)) > blockcutoff || sum(M(end-intrusion*blocksize+1:end,nn)) < sum(M(1:intrusion*blocksize,right))
    right = 0;
end



    if nn > Im_size(1)  && left ~= 0
       fprintf('%d on %d\n',nn,left)
    residual = [X(last_el-intrusion*blocksize+1:last_el,left);X(1:(blocksize-intrusion)*blocksize,nn)];
    Mloc = [M(last_el-intrusion*blocksize+1:last_el,left);M(1:(blocksize-intrusion)*blocksize,nn)];
    Mneighbour = M(last_el-intrusion*blocksize+1:last_el,left);
    Mneighbour = [ones(blocksize*(blocksize-intrusion),1);Mneighbour];
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
    residual = Mloc.*residual;
        toshow = X(:,left);
    %stencil left
        print = 10;
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
%         
%         if norm(residual)/print < 2
%         NeibourIntrudingPlotter([X(last_el-intrusion*blocksize+1:last_el,left);...
%             X(1:(blocksize-intrusion)*blocksize,nn)],U,Zloc,16)
%         print = norm(residual)/5;
%         end
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,left);
        FixedNeighbour = M(last_el-intrusion*blocksize+1:last_el,left);
        Xoriginal(Mneighbour==0) = Xloc(FixedNeighbour==0);
        
        X(:,left) = Xoriginal;
        M(last_el-intrusion*blocksize+1:last_el,left) = rand(intrusion*blocksize,1)>Mask_retention;    
    
%         NeighbourImprovement_check(toshow,X(:,left),16);pause;
%     if sum(M(:,left))/numel(M(:,left)) > masking_quality_cutoff seq(seq==left)=[];end
        
    end
    

    if nn <= Im_size(1)*(Im_size(2)-1) && right ~= 0
%     stencil right
    fprintf('%d on %d\n',nn,right)
    toshow = X(:,right);
    residual = [X(last_el-(blocksize-intrusion)*blocksize+1:last_el,nn);X(1:intrusion*blocksize,right)];
    Mloc = [M(last_el-(blocksize-intrusion)*blocksize+1:last_el,nn);M(1:intrusion*blocksize,right)];
    MofNeighbour = [M(1:intrusion*blocksize,right);ones(blocksize*(blocksize-intrusion),1)];
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
        
%         NeibourIntrudingPlotter([X(last_el-(blocksize-intrusion)*blocksize+1:last_el,nn);...
%             X(1:intrusion*blocksize,right)],U,Zloc,16)
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,right);
        GhostMask = [ones(blocksize*(blocksize-intrusion),1);M(1:intrusion*blocksize,right)];
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        X(:,right) = Xoriginal;
        M(1:intrusion*blocksize,right) = rand(intrusion*blocksize,1)>Mask_retention;
%         NeighbourImprovement_check(toshow,X(:,right),16);pause;
        
%     if sum(M(:,right))/numel(M(:,right)) > masking_quality_cutoff seq(seq==right)=[];end
    end

    
    
    

    if mod(nn,Im_size(1)) ~= 1 && up ~= 0
%     stencil up
    fprintf('%d on %d\n',nn,up)
    residual = zeros(length(X(:,nn)),1);
    residual(topintrusion) = X(bottomintrusion,up);
    residual(bottombase) = X(topbase,nn); 
    toshow = X(:,up);
    Mloc = zeros(length(M(:,nn)),1);
    Mloc(topintrusion) = M(bottomintrusion,up);
    Mloc(bottombase) = M(topbase,nn);
    MofNeighbour = ones(size(M,1),1);
    MofNeighbour(bottomintrusion) = M(bottomintrusion,up);
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
        
%         NeibourIntrudingPlotter(toshow,U,Zloc,16)
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,up);
        GhostMask = ones(last_el,1); %last element is the numel of the patch
        GhostMask(topintrusion) = Mloc(topintrusion);
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        X(:,up) = Xoriginal;
        M(bottomintrusion,up) = rand(intrusion*blocksize,1)>Mask_retention;
    
%         NeighbourImprovement_check(toshow,X(:,up),16);pause;
        
%     if sum(M(:,up))/numel(M(:,up)) > masking_quality_cutoff seq(seq==up)=[];end
    end 
    
    if mod(nn,Im_size(1)) ~= 0 && down ~= 0
%     stencil down
    fprintf('%d on %d\n',nn,down)
    
%     topsection = 1:blocksize;
%     topsection = topsection-1;
%     topsection = topsection * blocksize;
%     topsection = repmat(topsection,(intrusion),1);
%     topsection = topsection + repmat([1:intrusion]',1,blocksize);
%     topsection = reshape(topsection,numel(topsection),1);
% 
%     bottomsection = blocksize*blocksize - (repmat(fliplr(1:blocksize)-1,blocksize-intrusion,1)*blocksize + repmat(fliplr(0:blocksize-intrusion-1)',1,blocksize));

    residual = zeros(length(X(:,nn)),1);
    residual(topbase) = X(bottombase,nn);
    residual(bottomintrusion) = X(topintrusion,down); 
    toshow = X(:,down);
    Mloc = zeros(length(M(:,nn)),1);
    Mloc(topbase) = M(bottombase,nn);
    Mloc(bottomintrusion) = M(topintrusion,down);
    MofNeighbour = ones(size(M,1),1);
    MofNeighbour(topintrusion) = M(topintrusion,down);
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
%                 NeibourIntrudingPlotter(toshow,U,Zloc,16)
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,down);
        GhostMask = ones(last_el,1); %last element is the numel of the patch
        GhostMask(bottomintrusion) = Mloc(bottomintrusion);
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        X(:,down) = Xoriginal;
        M(topintrusion,down) = rand(intrusion*blocksize,1)>Mask_retention;
    
%         NeighbourImprovement_check(toshow,X(:,down),16);pause;
        
%     if sum(M(:,down))/numel(M(:,down)) > masking_quality_cutoff seq(seq==down)=[];end
    
    end
    
end

