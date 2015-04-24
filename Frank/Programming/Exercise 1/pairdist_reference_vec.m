function D = pairdist_reference_vec(P, Q)

[X1, X2] = meshgrid(P(:,1), Q(:,1));
[Y1, Y2] = meshgrid(P(:,2), Q(:,2));
D = sqrt((X1 - X2).^2 + (Y1 - Y2).^2)';
end
