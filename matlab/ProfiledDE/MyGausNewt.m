function newcoefs = MyGausNewt(coefs,lsopts,basis_cell,Ycell,Tcell,...
    loads,lambda,fn,pars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function MyGausNewt
%
% Gauss-Newton Estimation for use with SplineCoefErrcell.
% Not currently used or recommended. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[res0,Dres] = SplineCoefErrcell(coefs,basis_cell,Ycell,Tcell,loads,...
    lambda,fn.fn,fn.dfdx,pars);


coef0 = coefs;
F0 = sum(res0.^2);
F1 = 0;
iter = 0;
fundif = 1;
gradnorm0 = mean(res0'*Dres);
gradnorm1 = 1;

while gradnorm1 > lsopts.TolFun | fundif > lsopts.TolFun
    iter = iter + 1;
    if iter > lsopts.MaxIter, break; end
    Dcoef = Dres\res0;
    coef1 = coef0 - Dcoef;
    [res1,Dres] = SplineCoefErrcell(coef1,basis_cell,Ycell,Tcell,loads,...
        lambda,fn.fn,fn.dfdx,pars);
    F1 = sum(res1.^2);
    gradnorm1 = mean(res1'*Dres);
    fundif = abs(F0-F1)/abs(F0);
    if isequal(lsopts.Display,'iter')
        disp(num2str([iter, F1, fundif, gradnorm1]));
    end
    coef0 = coef1;
    res0 = res1;
    F0 = F1;
    gradnorm0 = gradnorm1;
end

newcoefs = coef0;

end