function [f,J] = JointSplineErr(bcoefs,ncoef,basis_cell,Ycell,Tcell,wts,lambda,fn,alg,fn_extras)

J = [];

coefs = bcoefs(1:ncoef);
pars = bcoefs((ncoef+1):length(bcoefs));


if nargout > 1,

    [f,J1] = SplineCoefErr(coefs,basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,fn_extras);


    DEfd = Make_fdcell(coefs,basis_cell);

    [f2,J2] = SplineParsErr(pars,DEfd,fn,fn_extras);
    
    J =  [J1 [zeros(size(J1,1)-size(J2,1),size(J2,2)); J2] ];
            
else

    f = SplineCoefErr(coefs,basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,fn_extras);

end


end
