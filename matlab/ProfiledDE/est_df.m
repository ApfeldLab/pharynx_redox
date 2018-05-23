function [dhYdY,dpdY,df] = est_df(DEfd,fn,Ycell,Tcell,wts,lambda,pars,alg,fn_extras)

if nargin<9, fn_extras = []; end
if nargin<8, alg = []; end

d2Gdc2 = make_d2gdc2(DEfd,fn,Tcell,wts,lambda,pars,alg,fn_extras);

d2Gdcdp = make_d2gdcdp(DEfd,fn,lambda,pars,alg,fn_extras);

[f,Zmat] = djdc_pts(DEfd,Ycell,Tcell,wts);

d2Jdp2 = make_d2jdp2(DEfd,fn,Ycell,Tcell,lambda,pars,alg,wts,fn_extras);

d2JdpdY = make_d2jdpdy(DEfd,fn,Ycell,Tcell,lambda,pars,alg,wts,fn_extras);

dpdY = -d2Jdp2\d2JdpdY;

dhYdY = Zmat*(d2Gdc2\(d2Gdcdp*dpdY + Zmat'));

I = eye(size(dhYdY,1));

df = [trace((I-dhYdY)'*(I-dhYdY)) trace(dhYdY) trace(dhYdY'*dhYdY)];