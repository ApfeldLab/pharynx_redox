function d2J = make_d2jdcdy(fd_cell,Tcell,wts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_d2jdcdy
%
% Calculates the second derivative of squared error with respect to
% coefficients and data. 
%
% INPUTS:
%
% fd_cell    - a cell array, each component contains a functional data 
%                object corresponding to the same component of the DE. 
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% if wts is empty, the inverse standard deviation of each component
% is used as a weighting for each dimension. If lambda is a singleton, it 
% is multiplied by the average weight in each component. 
%
% OUTPUT:
%
% d2J   - second derivative of squared error with respect to coefficients
%         and data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lambda = 0;
Ycell = eval_fdcell(Tcell,fd_cell,0);
[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);

basis_cell = getcellbasis(fd_cell);

Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);

for i = 1:length(Zmat_cell)
    if ~isempty(Tcell{i})
        n = length(wts(i));
        Zmat_cell{i} = spdiags(wts{i},0,n,n)*Zmat_cell{i};
    else
        Zmat_cell{i} = zeros(0,getnbasis(basis_cell{i}));
    end
end
d2J = -mattdiag_cell(Zmat_cell,0)';


