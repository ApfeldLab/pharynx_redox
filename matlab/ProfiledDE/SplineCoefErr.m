 function [f,J] = SplineCoefErr(coefs,basis_cell,Ycell,Tcell,wts,lambda,...
      fn,alg,pars,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SplineCoefErr
%
% error and jacobian calculations for fitting a spline to data with a DE
% penalty
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

if nargin<10
    fn_extras = [];
end


% Convert wts, lambda and Ycell to the right format:

[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

lambda = sqrt(lambda);

% By default, system has no algebraic components. 

if isempty(alg) 
    alg = ones(length(Ycell),1);
end

% Re-pass parameter-coefficient matrix:

newcoefs = coefs;

% Create the appropriate functional data object

fd_cell = Make_fdcell(newcoefs,basis_cell);


% Matrix of penalty at quadrature points

quadvals = getquadvals(basis_cell{1});
qvals = sqrt(quadvals(:,2));
m = length(qvals);

dpredns = cell2mat(eval_fdcell(quadvals(:,1),fd_cell,alg));

fvals = cell2mat(feval_check(fn.fn,quadvals(:,1),fd_cell,pars,fn_extras));
dfvals = feval_check(fn.dfdx,quadvals(:,1),fd_cell,pars,fn_extras);


pen = (dpredns - fvals);

pen = kron(qvals,lambda).*pen;

pen = reshape(pen,size(pen,1)*size(pen,2),1);


% Now to calculate a Jacobian

if nargout > 1
    
    [err,Zmat] = djdc_pts(fd_cell,Ycell,Tcell,wts);

    
    % Derviatives of penalty wrt coefs:

    Dphimat_cell = getvalues_cell(basis_cell,alg);
    phimat_cell = getvalues_cell(basis_cell,0);
    
    dpen2 = cell(length(phimat_cell));
    
    for i = 1:length(Ycell)
        Dphimat_cell{i} = spdiags(qvals,0,m,m)*...
            sparse(lambda(i)*Dphimat_cell{i});
        for j = 1:length(Ycell)
            dpen2{i,j}=spdiags(dfvals{i,j}.*qvals,0,m,m)*...
                sparse(lambda(i)*phimat_cell{j});
        end
        
    end

    dpen1 = mattdiag_cell(Dphimat_cell,0);
    dpen2 = cell2mat(dpen2);
    
    J1 = [-Zmat; (dpen1-dpen2)];
 
    J = J1;
   
else    
    
    err = djdc_pts(fd_cell,Ycell,Tcell,wts);
    
end

    % Putting the error vector together. 
    
    f = [err; pen];
        
end
