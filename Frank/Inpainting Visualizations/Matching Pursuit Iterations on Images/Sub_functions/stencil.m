%stencil

function [X,M,seq] =  stencil(X,M,U,nn,blocksize,Im_size,seq,masking_quality_cutoff,neibounorm,rc_min)
Mask_retention = .5; % 1 for all ones
lowbound = 0.0;
intrusion = floor(12);

last_el = blocksize^2;
blockcutoff = masking_quality_cutoff * intrusion * blocksize;


neibounorm = neibounorm *1;
rc_min = rc_min*1 ;


[upintrude,uphost,downintrude,downhost] = neighbourhood_layout_up_down(blocksize,intrusion);

%%
up = nn-1;
down = nn+1;
left = nn - Im_size(1);
right = nn + Im_size(1);

[up,down,left,right] = check_neighbour(seq,up,down,left,right,...
    blockcutoff,lowbound,blocksize,intrusion,upintrude,downintrude,M,nn);

%%
    if nn > Im_size(1)  && left ~= 0
%        fprintf('%d on %d\n',nn,left)
    residual = [X(last_el-intrusion*blocksize+1:last_el,left);X(1:(blocksize-intrusion)*blocksize,nn)];
    Mloc = [M(last_el-intrusion*blocksize+1:last_el,left);M(1:(blocksize-intrusion)*blocksize,nn)];
    Mneighbour = M(last_el-intrusion*blocksize+1:last_el,left);
    Mneighbour = [ones(blocksize*(blocksize-intrusion),1);Mneighbour];
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
    residual = Mloc.*residual;
    %stencil left
        it = 0;
        rc_max = 1;
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
        vec1 = (Uloc'*(residual));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
        it = it+1;
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,left);
        FixedNeighbour = M(last_el-intrusion*blocksize+1:last_el,left);
        Xoriginal(Mneighbour==0) = Xloc(FixedNeighbour==0);
        M(last_el-intrusion*blocksize+1:last_el,left) = 1-((1-M(last_el-intrusion*blocksize+1:last_el,left))...
        .*rand(intrusion*blocksize,1)>Mask_retention);    
        X(:,left) = Xoriginal.*M(:,left);
    if sum(M(:,left))/numel(M(:,left)) > masking_quality_cutoff seq(seq==left)=[];end
        
    end
    

    if nn <= Im_size(1)*(Im_size(2)-1) && right ~= 0
%     stencil right
%     fprintf('%d on %d\n',nn,right)
    toshow = X(:,right);
    residual = [X(last_el-(blocksize-intrusion)*blocksize+1:last_el,nn);X(1:intrusion*blocksize,right)];
    Mloc = [M(last_el-(blocksize-intrusion)*blocksize+1:last_el,nn);M(1:intrusion*blocksize,right)];
    MofNeighbour = [M(1:intrusion*blocksize,right);ones(blocksize*(blocksize-intrusion),1)];
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
    it = 0;
    rc_max = 1;
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)

        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
        
        it = it + 1;
        end 
                
        
        Xloc = U*Zloc;
        Xoriginal = X(:,right);
        GhostMask = [ones(blocksize*(blocksize-intrusion),1);M(1:intrusion*blocksize,right)];
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        
        M(1:intrusion*blocksize,right) = 1-((1-M(1:intrusion*blocksize,right))...
            .*rand(intrusion*blocksize,1)>Mask_retention);
        
        X(:,right) = Xoriginal.*M(:,right);
    

    if sum(M(:,right))/numel(M(:,right)) > masking_quality_cutoff seq(seq==right)=[];end
    end

    
    
    

    if mod(nn,Im_size(1)) ~= 1 && up ~= 0
%     stencil up

    residual = zeros(length(X(:,nn)),1);
    residual(upintrude) = X(downintrude,up);
    residual(downhost) = X(uphost,nn); 
    toshow = X(:,up);
    Mloc = zeros(length(M(:,nn)),1);
    Mloc(upintrude) = M(downintrude,up);
    Mloc(downhost) = M(uphost,nn);
    MofNeighbour = ones(size(M,1),1);
    MofNeighbour(downintrude) = M(downintrude,up);
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
        it = 0 ;
        rc_max = 1;
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)

        vec1 = (Uloc'*(residual));%./sqrt(sum((repmat(m,1,length(residual)).*U).^2,2));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);
        
%         NeibourIntrudingPlotter(toshow,U,Zloc,16)
        it = it + 1;
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,up);
        GhostMask = ones(last_el,1); %last element is the numel of the patch
        GhostMask(upintrude) = Mloc(upintrude);
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        
        M(downintrude,up) = 1-((1-M(downintrude,up))...
        .*rand(intrusion*blocksize,1)>Mask_retention);
        X(:,up) = Xoriginal.*M(:,up);

    if sum(M(:,up))/numel(M(:,up)) > masking_quality_cutoff seq(seq==up)=[];end
    end 
    
    if mod(nn,Im_size(1)) ~= 0 && down ~= 0
%     stencil down
    residual = zeros(length(X(:,nn)),1);
    residual(uphost) = X(downhost,nn);
    residual(downintrude) = X(upintrude,down); 
    toshow = X(:,down);
    Mloc = zeros(length(M(:,nn)),1);
    Mloc(uphost) = M(downhost,nn);
    Mloc(downintrude) = M(upintrude,down);
    MofNeighbour = ones(size(M,1),1);
    MofNeighbour(upintrude) = M(upintrude,down);
    Uloc = U.*(repmat(Mloc,1,size(U,2)));
    Zloc = zeros(size(Uloc,2),1);
        it = 0;
        rc_max = 1;
        while (norm(residual) >= neibounorm*norm(X(:,nn)) && rc_max > rc_min)
%         while (norm(residual) >= neibounorm*norm(X(:,nn)) || rc_max > rc_min) && it < it_stop
        vec1 = (Uloc'*(residual));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Zloc(arg) = Zloc(arg) + vec1(arg)/norm(d);

        it = it +1;
        end 
        Xloc = U*Zloc;
        Xoriginal = X(:,down);
        GhostMask = ones(last_el,1); %last element is the numel of the patch
        GhostMask(downintrude) = Mloc(downintrude);
        Xoriginal(MofNeighbour==0) = Xloc(GhostMask==0);
        
        M(upintrude,down) = 1-((1-M(upintrude,down))...
        .*rand(intrusion*blocksize,1)>Mask_retention);
        X(:,down) = Xoriginal.*M(:,down);
        
        
    if sum(M(:,down))/numel(M(:,down)) > masking_quality_cutoff seq(seq==down)=[];end
    
    end
  
    
end

