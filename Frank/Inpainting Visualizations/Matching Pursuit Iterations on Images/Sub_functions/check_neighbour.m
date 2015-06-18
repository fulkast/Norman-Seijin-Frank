function[up,down,left,right] = check_neighbour(seq,up,down,left,right,blockcutoff...
                                ,lowbound,blocksize,intrusion,upintrude,downintrude,...
                                M,nn)
               
if ~sum(up==seq)
    up = 0;
elseif sum(M(downintrude,up)) >= blockcutoff || sum(M(downintrude,up)) < lowbound || sum(M(upintrude,nn)) < sum(M(downintrude,up))
    up = 0;
end


if ~sum(down==seq)
    down = 0;
elseif sum(M(upintrude,down)) >= blockcutoff || sum(M(upintrude,down)) < lowbound || sum(M(downintrude,nn)) < sum(M(upintrude,down))
    down = 0;
    
end
 
if ~sum(left==seq) 
    left = 0; 
elseif sum(M(end-intrusion*blocksize+1:end,left)) >= blockcutoff  || sum(M(upintrude,left)) < lowbound || sum(M(1:intrusion*blocksize,nn)) < sum(M(end-intrusion*blocksize+1:end,left))
    left = 0;
end

if ~sum(right==seq) 
    right = 0;
elseif  sum(M(1:intrusion*blocksize,right)) >= blockcutoff  || sum(M(upintrude,right)) < lowbound || sum(M(end-intrusion*blocksize+1:end,nn)) < sum(M(1:intrusion*blocksize,right))
    right = 0;
end

end