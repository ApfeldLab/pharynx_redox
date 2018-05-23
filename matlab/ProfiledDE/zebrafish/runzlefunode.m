function y = heartfunode(t,r,p)

y = r;

y(1) = p(1) + p(2)*r(1) + p(3)*r(1).^2 + p(4)*r(1).^3 + p(5)*r(2) + p(6)*r(1).*r(2);
y(2) = p(7) + p(8)*r(1) + p(9)*r(2);
