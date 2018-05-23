function [coefs,resnorm,resid,eflag,output,J] = SplineCoef_lsq(coefs,lsopts_in,...
        basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,betadef)
    
eflag = 0;
iter = 1;
    
while eflag <= 0 && iter < lsopts_in.maxiter2
    [coefs,resnorm,resid,eflag,alpha,output,J] = lsqnonlin(@SplineCoefErr,coefs,[],[],lsopts_in,...
         basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,betadef);
    iter = iter+1;
end