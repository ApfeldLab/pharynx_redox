function [f,J] = SplineCoefErr_DEfit(coefs,DEfd,ind,fn,pars,alg,fn_extras)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function SplineParsErr
%
% Residuals and derivatives for fitting parameters to a DE smooth. 
%
% INPUTS
%
% coefs -   current estimated co-efficients for the components of interest
%
% DEfd  -   a cell array of functional data objects representing the smooth
%           of data to a DE.
%
% ind   -   a vector of indicies indicating which components are to be
%           re-fit to the differential equation, with the others fixed. 
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
% pars  -   current parameters to the differential equation
%
% fn_extras - object containing additional input to the right hand side 
%               of the differential equation. 
%
% OUTPUT:
%
% f  -  vector residuals of the fit, and of the penalty at quadrature
%       points.
%
% J  -  derivatives of f with respect to parameters.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<7
    fn_extras = [];
end

if isempty(alg) 
    alg = ones(length(path),1);
end

basis_cell = getcellbasis(DEfd);

DEfd = update_fdcell(coefs,ind,DEfd);

quadvals = getquadvals(basis_cell{1});
qvals = sqrt(quadvals(:,2));
m = length(qvals);

dpredns = cell2mat(eval_fdcell(quadvals(:,1),DEfd,1));

fvals = cell2mat(feval_check(fn.fn,quadvals(:,1),DEfd,pars,fn_extras));

pen = (dpredns - fvals);

pen = kron(qvals,ones(1,(size(pen,2)))).*pen;

f = reshape(pen,size(pen,1)*size(pen,2),1);

if nargout > 1

   dfvals = feval_check(fn.dfdx,quadvals(:,1),DEfd,pars,fn_extras);
   Dphimat_cell = getvalues_cell(basis_cell,alg);
   phimat_cell = getvalues_cell(basis_cell,0);
   

   J = cell(size(dfvals));
   J(:) = {0};
   
   for i = 1:length(DEfd)
        J{i,i} = J{i,i} + spdiags(qvals,0,m,m)*sparse(Dphimat_cell{i});
        for j = 1:length(DEfd)
           J{i,j} = J{i,j} - spdiags(dfvals{i,j}.*qvals,0,m,m)*...
                sparse(phimat_cell{j});
        end
        
   end

   J = cell2mat(J(:,ind));
    
end

end