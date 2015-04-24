function Y = stdize(X)
Y = (X - repmat(mean(X),[length(X) 1])) ./ repmat(std(X),[length(X) 1]);
end

