function [yout,dydx] = symmetrizedSymmetricDirichlet(J,triAreas)

a = J(1,1,:);
b = J(1,2,:);
c = J(2,1,:);
d = J(2,2,:);
detJ = a.*d-b.*c;
invJ = [d -b; -c a]./detJ;
a2 = invJ(1,1,:);
b2 = invJ(1,2,:);
c2 = invJ(2,1,:);
d2 = invJ(2,2,:);

y1 = (a.*a + b.*b + c.*c + d.*d);
y2 = (a2.*a2 + b2.*b2 + c2.*c2 + d2.*d2);

y = (y1 + y2).*(abs(detJ)+1);
y = y - 8; 
yout = dot(y(:),triAreas);
dydx = dlgradient(yout,J);

end