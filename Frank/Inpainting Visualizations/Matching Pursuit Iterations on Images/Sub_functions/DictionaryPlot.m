function imnew = DictionaryPlot(U,Im_size,d)

side = d;

mycell = mat2cell(U,d^2,ones(Im_size(1)*Im_size(2),1));
% nrow = round(sqrt(size(U,2)));

mycell = reshape(mycell,Im_size(1),Im_size(2));

backfun = @(mycell) reshape(mycell,side,side);
imnew = cell2mat(cellfun(backfun,mycell,'UniformOutput',false));


end