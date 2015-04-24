function [ q ] = extract( I,d )
%EXTRACT Summary of this function goes here
%   Detailed explanation goes here

[r, c ,v] = size(I);
q = zeros(d^2,ceil(r/d)*ceil(c/d),v);
padding = @(s) padarray(s,[d-size(s,1) d-size(s,2)],'replicate','post');
ourfun = @(imcell) reshape(imcell,d^2,1);


for i = 1: size(I,3)

imcell = mat2cell(I(:,:,i),[d*ones(1,floor(size(I(:,:,i),1)/d)) mod(size(I(:,:,i),1),d)],[d*ones(1,floor(size(I(:,:,i),2)/d)) mod(size(I(:,:,i),2),d)]);

% X = cellfun(@mean2,imcell);

if isempty(imcell{end,1})
    imcell = imcell(1:end-1,:);
end
   
if isempty(imcell{1,end})
    imcell = imcell(:,1:end-1);
end



imcell(:,end) = cellfun(padding,imcell(:,end),'UniformOutput',false);
imcell(end,:) = cellfun(padding,imcell(end,:),'UniformOutput',false);
X = cell2mat(cellfun(ourfun,imcell,'UniformOutput',false));
X = reshape(X,d^2,numel(imcell));

q(:,:,i) = X;
end

end

