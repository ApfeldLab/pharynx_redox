function val = getvalues_cell(basis_cell,k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function getvalues_cell
%
% Gets the stored values for a cell array of basis objects, basis_cell, 
% with derivative k. 
%
% k may be a singleton, or it may be an array, each element corresponding
% to one element of basis_cell. 
%
% Outputs a cell array of matrices each element corresponding to the values
% obtained from the corresponding element of basis_cell. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


val = cell(size(basis_cell));

if length(k)==1  
    k = k*ones(size(basis_cell));
end

for i = 1:numel(basis_cell)
    val{i} = getvalues(basis_cell{i},k(i));
end

end