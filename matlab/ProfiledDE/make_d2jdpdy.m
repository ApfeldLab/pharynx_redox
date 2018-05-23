function d2J = make_d2jdpdY(fd_cell,fn,Ycell,Tcell,lambda,pars,alg,wts,...
    fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdpdY
%
% The second derivative of squared error with respect to the parameters and
% observations.
%
% INPUTS
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% fn    -   a struct containing functions of the right hand side and its
%           derivatives. Should include at least:
%
%  fn.fn:       a handle to a function of the form f(t,y,p) giving the
%               right hand side of the DE, should accept vector t, cell
%               array of function data object, fd_cell, and output cell 
%               array vectors, one for each component. 
%
%  fn.dfdx:     a handle to a function giving the derivative of fn.fn with 
%               respect to y. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.d2fdx2:   a handle to a function giving the Hessian of fn.fn with
%               respect to y. Output should be as above, with derivatives 
%               given in the second dimension of the cell array. 
%
%  fn.d3fdx3:   a handle to a function giving the third derivative of fn.fn
%               2nd dimnesion) thrice with respect to y. 
%
%  fn.dfdp:     a handle to a function giving the derivative of fn.fn with 
%               respect to p. Output should be a cell array of vectors, the
%               first dimenion giving the components of fn.fn with
%               derivatives indexed by the second dimension. 
%
%  fn.d2fdxdp:  a handle to  a function giving the second derivative of 
%               fn.fn with respect to y and p. The first dimension sound
%               index components of fn.fn, with subsequent dimensions
%               indexing derivatives with respect to x and p in that order.
%
%  fn.d3fdx2dp: a handle to a function giving the third derivative of fn.fn
%               with respect to y twice (second and third dimensions) and p
%               (fourth dimension). 
%
% Ycell  -  a cell array of observations corresponding to the observation
%           times given in Tcell, below. 
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components.
%
% lambda  - a vector of smoothing parameters corresponding to each
%           component of the DE.
%
% pars    - vector of the current parameters in the penalty.
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation, but with
%           no middle terms. 
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% OUTPUT:
%
% d2J   -  a matrix giving the second derivative of squared error with
%          respect to the parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<9
    fn_extras = [];
end

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

% First make the matrices

d2Gdcdp = make_d2gdcdp(fd_cell,fn,lambda,pars,alg,fn_extras);
d2Gdc2 = make_d2gdc2(fd_cell,fn,Tcell,wts,lambda,pars,alg,fn_extras);

d3Gdc3 = make_d3gdc3(fd_cell,fn,lambda,pars,alg,fn_extras);

d3Gdc2dp = make_d3gdc2dp(fd_cell,fn,lambda,pars,alg,fn_extras);

dJdc = make_djdc(fd_cell,Ycell,Tcell,wts);
d2Jdc2 = make_d2jdc2(Tcell,wts,fd_cell);

d2JdcdY = make_d2jdcdy(fd_cell,Tcell,wts);


dcdp = -d2Gdc2\d2Gdcdp;
dcdY = -d2Gdc2\d2JdcdY;


% Now put it all together

d2cdpdY = zeros([size(dcdp) size(dcdY,2)]);
d2JdpdY = zeros([size(dcdp,2) size(dcdY,2)]);

for i = 1:length(d3Gdc3)
    d2cdpdY(i,:,:) = dcdp'*d3Gdc3{i}*dcdY;
end


for i = 1:length(pars)
    d2cdpdY(:,i,:) = squeeze(d2cdpdY(:,i,:)) + d3Gdc2dp{i}'*dcdY;
    d2JdpdY(i,:) = -dJdc'*(d2Gdc2\squeeze(d2cdpdY(:,i,:)));
end


d2J = 2*(dcdp'*(d2JdcdY+d2Jdc2*dcdY) + d2JdpdY); 

end