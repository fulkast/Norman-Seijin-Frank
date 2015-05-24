function D = pairdist(P, Q)
[X1, X2] = meshgrid(Q(:,1), P(:,1));
[Y1, Y2] = meshgrid(Q(:,2), P(:,2));
D = sqrt((X1-X2).^2 + (Y1-Y2).^2);
end
