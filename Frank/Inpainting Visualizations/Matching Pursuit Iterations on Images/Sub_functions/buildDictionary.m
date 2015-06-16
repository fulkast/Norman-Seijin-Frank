function U = buildDictionary(dim , DCT)

% Builds a dictionary with atoms of specified dimension
%
% INPUT
% dim: The dimensionality of the dictionary atoms
%
% OUTPUT:
% U (d x l) dictionary with unit norm atoms



try 
    if ~DCT
    temp = load('dictionary.mat');
    U = temp.U;
%     V1=overDCTdict(dim,floor(dim));
%     U = [U, V1];
    display('Custom Dictionary Utilized');
    else
    temp = load('nullptr.mat');
    U = temp.U;    
    end
    
catch
    % Input the alternative here
   V1=haarTrans(dim); outtext = 'Haar Used';
%     V1=overDCTdict(dim,floor(dim)); outtext = 'DCT Used' ;
    U=[V1];%,V1];
    fprintf('%s\n',outtext);
end

   