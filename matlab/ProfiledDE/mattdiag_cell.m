function vals = mattdiag_cell(matt_cell,fill)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function mattdiag_cell
%
% Creates a large array out of a cell array containing arrays.
%
% INPUTS:
%
% matt_cell - a cell vector of matrices
%
% fill  - should the elements of matt_cell be placed on the diagonal (0) of
%         a large array, or should they be replicated down columns (1)?
%
% OUTPUT
%
% vals - a large array containing the elements of matt_cell either in the
%        diagnonal elements, or repeated down columns. Some regularity in
%        the elments of matt_cell is assumed. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vals = matt_cell{1};

if fill 
    o1 = ones(1,length(matt_cell));
    vals = kron(o1,vals);
    for i = 2:length(matt_cell)
        vals = [vals; kron(o1,matt_cell{i})];
    end
else
    for i = 2:length(matt_cell)
        newvals = matt_cell{i};
        vals = [vals sparse(size(vals,1),size(newvals,2));...
            sparse(size(newvals,1),size(vals,2)) newvals];
    end
end

end