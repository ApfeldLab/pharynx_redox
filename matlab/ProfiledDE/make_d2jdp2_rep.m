function d2J = make_d2jdp2_rep(fd_cell,fn,Ycell,Tcell,lambda,pars,...
    parind,alg,wts,fn_extras,d2pen,pen_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdp2_rep
%
% The second derivative of squared error with respect to the parameters
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
%  fn.d2fdp2:   a handle to a function giving the second derivative of
%               fn.fn respect to p. It should output a cell array as above.
%               
%  fn.d3fdx3:   a handle to a function giving the third derivative of fn.fn
%               thrice with respect to y. 
%
%  fn.d3fdx2dp: a handle to a function giving the third derivative of fn.fn
%               with respect to y twice (second and third dimensions) and p
%               (fourth dimension). 
%
%  fn.d3fdxdp2: a handle to a function giving the third derivative of fn.fn
%               with respect to y (second dimension) and p twice (third and
%               fourth dimensions). 
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
% parind  - A matrix; rows provide the indices of the parameters to be used
%           in each replication. 
%
% alg     - a vector to describe which variables are governed by
%           differential (entry = 1), as opposed to algebraic equations 
%           (entry = 0). All ones by default. Higher-valued integer terms
%           will represent a correspondingly high-order equation.
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% d2pen   - a function giving the second derivative of a penalty on the 
%           parameters it should accept the parameters and pen_extras as
%           arguments
%
% pen_extras - additional input into d2pen
% 
% OUTPUT:
%
% d2J   -  a matrix giving the second derivative of squared error with
%          respect to the parameters. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<12
    pen_extras = [];
end

if nargin<11
    d2pen = [];
end

if nargin<10
    fn_extras = cell(size(fd_cell,1),1);
end

if isempty(parind)
    parind = repmat(1:length(pars),size(Ycell,1));
end

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

d2J = zeros(length(pars),length(pars));

for i = 1:size(fd_cell,1)
    d2Jdp2 = make_d2jdp2(fd_cell(i,:),fn,Ycell(i,:),Tcell(i,:),lambda(i,:),...
        pars(parind(i,:)),alg,wts(i,:),fn_extras{i});

    d2J(parind(i,:),parind(i,:)) = d2J(parind(i,:),parind(i,:)) + d2Jdp2;
end

if ~isempty(d2pen)
    d2J = d2J + d2pen(pars,pen_extras);
end

end