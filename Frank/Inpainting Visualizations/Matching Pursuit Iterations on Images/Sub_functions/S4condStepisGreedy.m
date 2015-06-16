% function to find the best atom
function [step1,innprod,step2,secondinnprod,conv,rc_max] = S4condStepisGreedy(U,patch_residual,conv)
% All inputs should be sent in after consideration of the masking

        
        t = cputime;
        vec1 = (U'*(patch_residual));
        vec = abs(vec1);
        Natoms = length(vec);
        conv = 0;
        [rc_max,traditional] = max(vec.*isfinite(vec));
        if (conv)
            step1 = traditional;
            innprod = vec1(step1(1));
            step2 = step1;
            secondinnprod = innprod;
        else
%         fprintf('naive move %.1d\n',traditional)
        residual = repmat(patch_residual,1,Natoms);
        
        projections = bsxfun(@times,U',vec1)';
        residual = residual - projections;
        
        vec2 = (U'*residual);
        
        %%longversion
        [~, secondmove] =  max(abs(vec2),[],1);
        ind = sub2ind(size(vec2),secondmove,1:Natoms);
        residual = residual - bsxfun(@times,U(:,secondmove),vec2(ind)) ;
        [proposedres ,step1 ] = min(sqrt(sum(residual.^2,1)));
        step2 = secondmove(step1);
%         cputime-t;
%         later... put in cell array and find the move combination with the
%         lowest residue
        innprod = vec1(step1(1));
        secondinnprod = vec2(step2(1),step1(1));
%         fprintf('Smart move %.1d\n',step1)
%         fprintf('Apparently the error will be %.1d\n',proposedres)
        conv = (step1 == traditional);
        %%longversionend
        
        %%quickversion
%         [step2, step1] = find(abs(vec2)==max(max(abs(vec2))));
%         step2 = step2(1);
%         step1 = step1(1);
%         innprod = vec1(step1);
%         secondinnprod = vec2(step2,step1);
%         
        %%quickversion
        
        
% pause
        end
  end