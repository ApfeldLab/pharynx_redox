function [f,J] = Spline_Lin(coefs,basiscell,Ycell,Tcell,wts,lambda,...
    fn,alg,pars,fn_extras)

[f,D] = SplineCoefErr(coefs,basiscell,Ycell,Tcell,wts,lambda,...
    fn,alg,pars,fn_extras);

f = f'*f;

J = f'*D;