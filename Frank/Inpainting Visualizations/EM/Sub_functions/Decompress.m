function I_rec = Decompress(I_comp)

% Your decompression code goes here!





[~, ~ ,v] = size(I_comp.U);
d = I_comp.d;
X = zeros(d^2,size(I_comp.Z,2),v);
for i = 1:v
X(:,:,i) = (I_comp.U(:,:,i)*I_comp.Z(:,:,i))+repmat(I_comp.mu(:,i),[1 size(I_comp.Z,2)]);

end


imnew = zeros(d*ceil(I_comp.aspect(1)/d),...
    d*ceil(I_comp.aspect(2)/d),v);

for i = 1:v
mycell = mat2cell(X(:,:,i),d^2,ones(1,size(X(:,:,i),2)));
mycell = reshape(mycell,ceil(I_comp.aspect(1)/d),...
    ceil(I_comp.aspect(2)/d));
backfun = @(mycell) reshape(mycell,d,d);
imnew(:,:,i) = cell2mat(cellfun(backfun,mycell,'UniformOutput',false));

end

imnew=imnew(1:I_comp.aspect(1),1:I_comp.aspect(2),:);
% figure
% imshow(imnew)

I_rec = imnew; % this is just a stump to make the evaluation script run, replace it with your code!