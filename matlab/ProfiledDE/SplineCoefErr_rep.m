function [f,J] = SplineCoefErr_rep(coefs,basis_cell,Ycell,Tcell,...
    wts,lambda,fn,alg,pars,parind,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SplineCoefErrcell_rep
%
% error and jacobian calculations for fitting a sequence of splines with a
% DE penalty. 
% 
% INPUTS:
% 
% coefs   -  co-efficients of the basis expansion
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
% parind -  an array of indices. The rows of parind provide the parameters
%           that apply to each replication of the experiment. 
%
% fn_extras   - object containing additional input to the right hand side 
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
    fn_extras = cell(size(basis_cell,1),1);
end

if isempty(parind)
    parind = repmat(1:length(pars),size(Ycell,1));
end

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

f = cell(size(basis_cell,1),1);
J = cell(size(basis_cell,1),1);

n1 = 1;
n2 = 1;

for i = 1:size(basis_cell,1)
    for j = 1:size(basis_cell,2)
        n2 = n2 + getnbasis(basis_cell{i,j});
    end

    if nargout > 1
        
        [f{i},J{i}] = SplineCoefErr(coefs(n1:(n2-1)),basis_cell(i,:),...
            Ycell(i,:),Tcell(i,:),wts(i,:),lambda(i,:),fn,alg,...
            pars(parind(i,:)),fn_extras{i});

    else
        f{i} = SplineCoefErr(coefs(n1:(n2-1)),basis_cell(i,:),...
            Ycell(i,:),Tcell(i,:),wts(i,:),lambda(i,:),fn,...
            alg,pars(parind(i,:)),fn_extras{i});
    end

    n1 = n2;
end

f = cell2mat(f);

if nargout>1
    J = mattdiag_cell(J,0);
end

end
