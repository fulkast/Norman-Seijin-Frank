% function to find the best atom
function [index,innprod,secondind,secondinnprod] = a7omfinder_4dynamicProgramming(U,patch_residual)
% All inputs should be sent in after consideration of the masking
        t = cputime;
        vec1 = (U'*(patch_residual));
        vec = abs(vec1);
        Natoms = length(vec);
        
        [~,traditional] = max(vec.*isfinite(vec));
        fprintf('naive move %.1d\n',traditional)
        residual = repmat(patch_residual,1,Natoms);
        projections = bsxfun(@times,U',vec1)';
        residual = residual - projections;
        
        vec2 = (U'*residual);
        residual2 = reshape(residual,numel(residual),1);
        residualcell2 = mat2cell(repmat(residual2,1,Natoms),Natoms*ones(Natoms,1),ones(Natoms,1));
        Us2s = mat2cell(repmat(U,Natoms,1),Natoms*ones(Natoms,1),ones(Natoms,1));
        veccell2 = num2cell(vec2);
        Us2s = cellfun(@times,Us2s,veccell2,'UniformOutput',false);
        residualcell2 = cellfun(@minus,residualcell2,Us2s,'UniformOutput',false);
        cputime-t
        
        % later... put in cell array and find the move combination with the
        % lowest residue
        
        fprintf('Smart move %.1d\n',index)
  end