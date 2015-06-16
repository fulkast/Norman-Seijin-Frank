function showDico(U,dim,Usize)
if (numel(find(Usize~=0))==0)
    Usize=[dim dim]
end
figure
imshow(0.5+my_col2im(0.5*(U-repmat((mean(U)),dim,1))./repmat(sqrt(var(U)),dim,1),sqrt(dim),Usize,0,[0 0]));

end