function [f,J] = djdc_pts(fd_cell,Ycell,Tcell,wts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function djdc_pts
%
% Calculates the point-wise errors between a cell array of functional data
% objects
%
% INPUTS:
%
% fd_cell    - a cell array, each component contains a functional data 
%                object corresponding to the same component of the DE. 
%
% Ycell  -  a cell array of observations, one element corresponding to one
%           measured component. 
%
% Tcell  -  vector of observation times. May be a vector or cell
%           array if times are different for different components. 
%
% wts     - the weighting for each observation of the DE. May be empty. If
%           a vector, assumed to be weights for each dimension. If a cell
%           array, gives a weighting per observation per dimension. 
%
% if wts is empty, all weights are assumed to be equal. 
%
% OUTPUT:
%
% f - vector of residuals of the functional data object from the Ycell
%
% J - derivative of these residuals with respect to co-efficients in
%     fd_cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lambda = 0;
[wts,lambda,Ycell] = weightslambdasYcell(wts,lambda,Ycell);


basis_cell = getcellbasis(fd_cell);

eYcell = eval_fdcell(Tcell,fd_cell,0);

err = cell(size(eYcell));

for i = 1:length(Tcell)
    if ~isempty(Tcell{i})
        err{i} = sqrt(wts{i}).*(Ycell{i} - eYcell{i});
    end
end

f = cell2mat(err');

if nargout > 1
       
    Zmat_cell = eval_basis_cell(Tcell,basis_cell,0);
    
    for i = 1:length(Zmat_cell)
        n = length(wts(i));
        if ~isempty(Tcell{i})
            Zmat_cell{i} = spdiags(sqrt(wts{i}),0,n,n)*Zmat_cell{i};
        else
            Zmat_cell{i} = zeros(0,getnbasis(basis_cell{i}));
        end
    end
    J = mattdiag_cell(Zmat_cell,0);
    

end