function d2J = make_d2jdc2(Tcell,wts,fd_cell);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdc2
%
% The second derivative of squared error with respect to the coefficients
% in the cell array of functional data objects fd_cell. 
%
% INPUTS
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components.
%
% wts   - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% fd_cell  -  a cell array of functional data objects giving the current
%               spline smooth to the data.
%
% OUTPUT:
%
% d2J   -   the second derivative matrix of the spline criterion with
%           respect to coefficients. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lambda = 0;
Ycell = eval_fdcell(Tcell,fd_cell,0);
[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

basis_cell = getcellbasis(fd_cell);

Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);

H = cell(length(Zmat_cell),1);
for i=1:length(basis_cell)
    if ~isempty(Tcell{i})
        m = length(wts{i});
        H{i} = Zmat_cell{i}'*spdiags(wts{i},0,m,m)*Zmat_cell{i};
    else
        H{i} = sparse(getnbasis(basis_cell{i}),getnbasis(basis_cell{i}));
    end
end
d2J = mattdiag_cell(H,0);


end