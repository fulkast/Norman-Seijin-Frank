function U = buildDictionary(dim)

% Builds a dictionary with atoms of specified dimension
%
% INPUT
% dim: The dimensionality of the dictionary atoms
%
% OUTPUT:
% U (d x l) dictionary with unit norm atoms



try 
    temp = load('dictionary.mat');
    U = temp.U;
catch
    % Input the alternative here
    V1=haarTrans(dim);
    V2=overDCTdict(dim,floor(dim));
    U=[V2];%,V1];
end

   