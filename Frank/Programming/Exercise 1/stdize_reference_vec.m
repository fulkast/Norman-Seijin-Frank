function Y = stdize_reference_vec(X)

m = mean(X);
s = std(X);
Y = (X - repmat(m,size(X,1),1)) ./ repmat(s,size(X,1),1);
end

