function [f,J,DEfd,dcdp] = ProfileErr(pars,DEfd,fn,lambda,Ycell,Tcell,...
    wts,alg,lsopts,fn_extras,pen,dpen,pen_extras,dcdp,oldpars)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ProfileErr
%
% A function to calculate a least squares error and gradient for profile
% fitting a differential equation
%
% INPUTS:
%
% pars  -   estimate of the parameters
%
% DEfd    - a cell array, each component contains a functional data 
%           object corresponding to the same component of the DE. This is
%           not entered, but assumed to be defined as a global variable.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:       a handle to a function of the form f(t,y,p) giving the 
%               right hand side of the DE, should accept vector t, a cell 
%               array of function data objects, fd_cell, and output cell 
%               cell array of vectors, one for each component. 
%
%  fn.dfdx:     a handle to a function giving the derivative of fn.fn with 
%               respect to y. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.dfdp:     a handle to a function giving the derivative of fn.fn with 
%               respect to p. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.d2fdx2:   a handle to a function giving the Hessian of fn.fn with
%               respect to y. Output should be as above, with derivatives 
%               given in the second dimension of the cell array. 
% 
%  fn.d2fdxdp:  a handle to  a function giving the second derivative of 
%               fn.fn with respect to y and p. The first dimension sound
%               index components of fn.fn, with subsequent dimensions
%               indexing derivatives with respect to x and p in that order.
%
% lambda  - a vector of smoothing parameters corresponding to each
%           component of the DE.
%
% Ycell    - the Ycell of observed values of the DE, either matrix or cell
%           array.
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% lsopts  - options for the profile least squares routine, may be empty. 
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% pen   - a penalty for the parameters; this should be a function that 
%         takes the parameters as arguments and additional input pen_extra
%
% dpen  - the derivative of the penalty with respect to parameters, this
%         should be a function that takes the parameters as arguments and
%         an additional input pen_extra
%
% pen_extras - additional input into penalty functions
%
% dcdp - derivative of co-efficients with respect to those parameters used
%        to estimate the current smooth. This is used to speed-up the
%        inner optimisation - no speed-up performed if empty or not
%        included. 
%
% oldpars - parameters used to estimate the current smooth. Useful for
%           speeding up the inner optimisation. No speedup is performed if
%           empty or not included. 
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting for each dimension. If lambda is a singleton, it 
% is multiplied by the average weight in each component. 
%
% OUTPUT:
%
% f - vector of residuals of the functional data object from the Ycell
%
% J - derivative of these residuals with respect to pars
%
% DEfd - the smooth at the current parameter values
%
% dcdp - the derivative of the smooth with respect to the parameters at the
%        current parameter values
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp(pars')

if nargin<15
    dcdp = [];
    oldpars = [];
end

if nargin<13
    pen_extras = [];
end

if nargin<12
    dpen = [];
end

if nargin<11
    pen = [];
end

if nargin<10
    fn_extras = [];
end

% Convert wts, lambda and Ycell to the right format:

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

    
% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(Ycell),1);
end

% Now lets start. 

coefs = getcellcoefs(DEfd);
basis_cell = getcellbasis(DEfd);


% Profile the co-efficients. 

if ~isempty(dcdp)
   coefs = coefs + dcdp*(pars-oldpars); 
end

if isequal(fn.fn,@genlinfun)
    newcoefs = genlin_smooth(Ycell,Tcell,wts,basis_cell,lambda,pars,alg,...
        fn_extras);
else
    newcoefs = lsqnonlin(@SplineCoefErr,coefs,[],[],lsopts,...
        basis_cell,Ycell,Tcell,wts,lambda,fn,alg,pars,fn_extras);
end

DEfd = Make_fdcell(newcoefs,basis_cell);


if nargout > 1 
       
    
    [f,Zmat] = djdc_pts(DEfd,Ycell,Tcell,wts);

    d2Gdc2 = make_d2gdc2(DEfd,fn,Tcell,wts,lambda,pars,alg,fn_extras);

    d2Gdcdp = make_d2gdcdp(DEfd,fn,lambda,pars,alg,fn_extras);
    
    dcdp = -d2Gdc2\d2Gdcdp;
    
    J = -Zmat*dcdp;
    
    if ~isempty(pen)      
        f = [f; pen(pars,pen_extras)];    
        
        J = [J; dpen(pars,pen_extras)];
    end
    
else
    
    f = djdc_pts(DEfd,Ycell,Tcell,wts);

    if ~isempty(pen)
        f = [f; pen(pars,pen_extras)];
    end
end

end
