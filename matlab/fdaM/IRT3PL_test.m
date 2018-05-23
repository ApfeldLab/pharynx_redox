n = 101;

theta = linspace(-3,3,n)';

Wbasis = create_IRT3PL_basis([-3,3],3);

plot(Wbasis)

Wmat = eval_basis(theta,Wbasis);

plot(theta, Wmat, '-')

coef = [0, 1, 1]';

Wfd = fd(coef, Wbasis);

plot(Wfd)

plot(Wfd,1)

plot(Wfd,2)

for nderiv=0:2
    Wpen = eval_penalty(Wbasis, nderiv);
    Wpen
    pause
end