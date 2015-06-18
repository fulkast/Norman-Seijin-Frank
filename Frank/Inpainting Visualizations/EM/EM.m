% EM code
function Z2 = EM(m,residual,Yobsloc,U,converged,sigma,rc_min,l)

% residual = Yobsloc +(1-m).*(U*Z(:,nn));

  while ~converged
    % Expectation  
    InitialResidual = residual;
    rc_max = max(abs(U'*residual));
    Z2 = zeros(l,1);
    % Maximization step
    while (norm(residual) > 0.1*sigma*norm(InitialResidual)) && (rc_max > rc_min)
        vec1 = (U'*(residual));
        vec = abs(vec1);
        [rc_max,arg] = max(vec.*isfinite(vec));
        d = U(:,arg);
        residual = residual - vec1(arg)/norm(d)*d;
        Z2(arg)=Z2(arg)+vec1(arg)/norm(d);
    end
    
    if norm((Yobsloc+(1-m).*(U*Z2))-InitialResidual)/norm(InitialResidual) < sigma
            converged = 1;
    end
        
    
    
  end