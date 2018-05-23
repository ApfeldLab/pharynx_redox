function basis_cell = getcellbasis(fd_cell)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getcellbasis
%
% Returns a cell array of basis objects corresponding to the input cell
% array of functional data objects fd_cell. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
basis_cell = cell(size(fd_cell));

for i = 1:numel(fd_cell) 
    basis_cell{i} = getbasis(fd_cell{i});
end

end