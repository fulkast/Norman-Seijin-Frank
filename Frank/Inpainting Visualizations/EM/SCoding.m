% SparseCoding.m
function Z = SCoding(Uloc,Z,residual,sigma,rc_min,nn)

InitialResidual = residual;
rc_max = max(abs(Uloc'*residual));

while (norm(residual) > norm(1*sigma*InitialResidual)) && (rc_max > 1*rc_min)
        vec1 = (Uloc'*(residual));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = Uloc(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Z(arg,nn)=Z(arg,nn)+vec1(arg)/norm(d);
        
end
    
