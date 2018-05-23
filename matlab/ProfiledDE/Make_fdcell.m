function fd_cell = Make_fdcell(coef,basis_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Make_fdcell
%
% A function to create a cell array of functional data objects.
%
% INPUTS:
%
% coef - a vector giving the co-efficients of the basis expansion, these
%        are concatenated across bases. 
% 
% basis_cell - a cell array of basis objects for coef to be applied to.
%
% OUTPUT:
%
% fd_cell - a cell array of functional data objects.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basis_cell = basis_cell';

fd_cell = cell(size(basis_cell));

here = 0;

for i = 1:numel(basis_cell)
    nbasis = getnbasis(basis_cell{i});
    fd_cell{i} = fd(coef((here+1):(here+nbasis)),basis_cell{i});
    here = here+nbasis;
end

fd_cell = fd_cell';

end