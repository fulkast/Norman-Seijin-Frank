function U = buildDictionary(dim)

% Builds a dictionary with atoms of specified dimension
%
% INPUT
% dim: The dimensionality of the dictionary atoms
%
% OUTPUT:
% U (d x l) dictionary with unit norm atoms


try 
    temp = load('projectedKmeanHS16pics.mat');
    U = temp.U;
catch
    % Input the alternative here
    warning('Building new default dictionary. Is that really okay???')
    V=overDCTdict(dim,floor(dim));
    U=[V];
end

   
