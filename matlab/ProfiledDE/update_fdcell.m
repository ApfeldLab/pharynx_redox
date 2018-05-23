function fd_cell = update_fdcell(coef,ind,fd_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Make_fdcell
%
% A function updates the co-efficients to some elements of a cell array of
% functional data obets. 
%
% INPUTS:
%
% coef - a vector giving the co-efficients of the basis expansion to be
%        updated. These are concatenated accross bases. 
%
% ind  - a vector of indeces indicating which co-efficients should be
%        updated. 
% 
% fd_cell - a cell array of functional data objects to be updated. 
%
% OUTPUT:
%
% fd_cell - a cell array of functional data objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis_cell = getcellbasis(fd_cell);

where = 0;
for i = 1:length(ind)
    fd_cell{ind(i)} = fd(coef(where + (1:getnbasis(basis_cell{ind(i)}))),...
        basis_cell{ind(i)});
    where = where + getnbasis(basis_cell{ind(i)});
end
