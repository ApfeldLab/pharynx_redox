function newpars = lintranspars(pars,s,k,a)

npars = pars;

npars = npars.*[s 1 1/s 1/s^2 s 1 1 1/s 1]';

newpars = npars;

newpars(1) = npars(1) - npars(2)*k + npars(3)*k^2 - npars(4)*k^3;
newpars(2) = npars(2) - 2*npars(3)*k + 3*npars(4)*k^2;
newpars(3) = npars(3) -3*npars(4)*k;
newpars(5) = npars(5) - npars(6)*k;
newpars(7) = npars(7) - npars(8)*k;

newpars = a*newpars;