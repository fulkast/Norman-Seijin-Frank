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
    L = 1000;
    U = [haarTrans(dim), overDCTdict(dim, L)];
end

   