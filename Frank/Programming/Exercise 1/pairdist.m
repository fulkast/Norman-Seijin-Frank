function D = pairdist(P, Q)
[thisx, thatx] = meshgrid(P(:,1),Q(:,1));
[thisy, thaty] = meshgrid(P(:,2),Q(:,2));

D = (sqrt((thisx-thatx).^2+(thisy-thaty).^2))';
end
