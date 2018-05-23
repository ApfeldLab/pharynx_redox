function newcoef = Spline_GausNewt(coefs,lsopts,basis_cell,Ycell,Tcell,wts,lambda,...
      fn,alg,pars,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Spline_GausNewt
%
% Gauss-Newton Minimization of SplineCoefErr
% 
% INPUTS:
% 
% coefs   -  co-efficients of the basis expansion
% 
% lsopts   - optimization options
% 
% basis_cell - a cell-array, each component should contain a basis for
%                the corresponding dimension of the DE. Assumes quadrature
%                values have been attached to each. 
% 
% Ycell   - observed values of the DE (matrix or cell array).
% 
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% lambda - vector of smoothing parameters for each dimension
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:   a handle to a function of the form f(t,y,p) giving the right 
%           hand side of the DE, should accept vector t, cell array 
%           function data object, fd_cell, and output cell array of
%           vectors, one for each component. 
%
%  fn.dfdx: a handle to a function giving the derivative of fn.fn with 
%           respect to y. Output should be a cell array of vectors, the
%           first dimenion giving the components of fn.fn with derivatives
%           indexed by the second dimension. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% pars  -   vector of inputs p into fn.fn. 
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting. If lambda is a singleton, it is multiplied by
% the average weight in each dimension. 
%
% OUTPUT:
%
% f  -  vector residuals of the fit, and of the penalty at quadrature
%       points.
%
% J  -  derivatives of f with respect to coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<11
    fn_extras = [];
end

% Convert wts, lambda and Ycell to the right format:

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

lambda = sqrt(lambda);

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(Ycell),1);
end

[res0,Dres] = SplineCoefErr(coefs,basis_cell,Ycell,Tcell,wts,lambda,...
           fn,alg,pars,fn_extras);
       
F0 = res0'*res0;
F1 = F0;

iter = 0;
fundif = 1;
gradnorm = 1;

H = speye(length(coefs));
g = zeros(length(coefs),1);

if isequal(lsopts.Display,'iter')
       header = sprintf(...
           '\n Iteration    Residual   Improvement    Grad-norm   ntry');
       disp(header)
       formatstr = ' %5.0f   %13.6g %13.6g %12.3g    %5.0f';
end

while gradnorm > lsopts.TolFun && fundif > lsopts.TolFun 
    
    iter = iter + 1;
    if iter > lsopts.MaxIter, break; end

 
    gnew = 2*Dres'*res0;
    
    gradnorm = norm(gnew);
    
    s = -H*gnew;
     
    ss = s/norm(s);
    
    %    alpha = 2;
%    trys = 0;
    disp([F0,F1])
    
%     alpha = lsqnonlin(@Spline_Lin,0,[],[],lsopts,coefs,ss,basis_cell,Ycell,Tcell,wts,lambda,...
%     fn,alg,pars,fn_extras);

    [alpha,fs,iter,errcode] = stepbhhh(@Spline_Lin,coefs,F0,gnew,s,25,basis_cell,Ycell,Tcell,wts,lambda,...
     fn,alg,pars,fn_extras);

    
%     while F0-F1 <= 0 && trys < lsopts.maxtry
%         alpha = 0.5*alpha;
%         tcoef = coefs + alpha*s;
%         [res0,Dres] = SplineCoefErr(tcoef,basis_cell,Ycell,Tcell,wts,lambda,...
%            fn,alg,pars,fn_extras);
%         F1 = res0'*res0;
%         trys = trys + 1;
%         fundif = F0 - F1;
%     end
    
    alpha = alpha/norm(s);
    coefnew = coefs + alpha*s;
    
  [res0,Dres] = SplineCoefErr(coefnew,basis_cell,Ycell,Tcell,wts,lambda,...
           fn,alg,pars,fn_extras);
       
    F1 = res0'*res0;
  
    
    talpha = alpha*norm(s)/norm(gnew); % Restart if H goes bad
    tcoef = coefs + talpha*gnew;
    tres = SplineCoefErr(tcoef,basis_cell,Ycell,Tcell,wts,lambda,...
           fn,alg,pars,fn_extras);

       

       
    y = (gnew-g)/alpha;
    s= alpha*s;
%    H = H + y*y'/(y'*s) - (H*s)*(H*s)'/(s'*H*s);

    H = H + (s*s')*(s'*y + y'*H*y)/(s'*y)^2 - (H*y*s' + s*(H*y)')/(s'*y);

       
    if tres'*tres < F1
        coefnew = tcoef;
        H = eye(length(gnew));
        F1 = tres'*tres;
        fundif = F0 - F1;
        trys = -1;
    end
    

    if isequal(lsopts.Display,'iter')
        iterstr = sprintf(formatstr,iter,F1,fundif,gradnorm,trys);
        disp(iterstr);
    end
    coefs = coefnew;
    g = gnew;
    F0 = F1;
end

newcoef = coefs;

end