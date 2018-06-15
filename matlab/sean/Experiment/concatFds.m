% Concatenate a bunch of fd objects together, assuming they have the same
% basis, etc.
% fdCells is an n x 1 cell array, containing a single fd object in each
% cell
function catFd = concatFds(fdCells)
    common_basis = getbasis(fdCells{1});
    coef_cells = cellfun(@(x) getcoef(x), fdCells, 'UniformOutput', false);
    
     % Check that coefficients are the same size
    coef_sizes = cellfun(@(x) size(x), coef_cells, 'UniformOutput', false);
    bases = cellfun(@(x) getbasis(x), fdCells, 'UniformOutput', false);
    
    all_coef_same = all(cellfun(@(x) isequal(coef_sizes{1}, x), coef_sizes));
    all_bases_same = all(cellfun(@(x) isequal(bases{1}, x), bases));
    
    if ~all_coef_same
        error('All FD coefficients must share the same dimensions');
    end
    if ~all_bases_same
        error('All FD must share the same basis');
    end
    
    coefs = cat(2, coef_cells{:});
    catFd = fd(coefs, common_basis);
end