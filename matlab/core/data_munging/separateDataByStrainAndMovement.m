% Suppress array growing warning, these cell arrays are tiny, and it makes
% the code much more readable
%#ok<*AGROW>


function [dataSeparated, labels] = separateDataByStrainAndMovement(data, indexer, region)
%SEPARATEDATABYINDEXERCOLUMNS Create cell array of data matrices for each
%   combination of (strain x movement), where movement falls in two
%   classes: (0 or 1) and (2 and 3).
%
%   The ordering of this array is as follows:
%       strain1, 0/1 movement
%       strain1, 2/3 movement
%       strain2, 0/1 movement
%       strain2, 2/3 movement
%       ...
%
%   This function also returns a cell array of labels (for plotting, etc.)
%   following the same order of the data cell array.
    
    strains = unique(indexer.Strain);
    dataSeparated = {};
    labels = {};
    
    low_mvmt_idx = indexer.(region) <= 1;
    hi_mvmt_idx  = indexer.(region) >= 2;
    
    for i=1:length(strains)
        strain = strains(i);
        
        strain_idx = indexer.Strain == strain;
        
        low_mvmt_strain_idx = low_mvmt_idx & strain_idx;
        hi_mvmt_strain_idx  = hi_mvmt_idx  & strain_idx;
        
        dataSeparated{end+1} = data(:, unique(indexer.Animal(low_mvmt_strain_idx)));
        dataSeparated{end+1} = data(:, unique(indexer.Animal(hi_mvmt_strain_idx)));
        
        labels{end+1} = strcat(strain, " (0/1)"); 
        labels{end+1} = strcat(strain, " (2/3)");
    end
end

