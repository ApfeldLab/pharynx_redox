function T = regionMeans(data)
%REGIONMEANS Create a table, with headers for each region, and each row
% corresponding to the mean of the values of a single animal inside that
% region.
%
% ARGUMENTS
%   data
%       An m x n matrix, where m = the size of a single profile measurement
%       and n = the number of animals (each column represents a single
%       animal).

    region_names = fieldnames(Constants.regions);
    n_worms = size(data, 2);
    n_regions = numel(region_names);
    
    varTypes = cell(1, n_regions);
    varTypes(:) = {'double'};
    T = table('Size', [n_worms n_regions], 'VariableTypes', varTypes, 'VariableNames', region_names);
    
    for i = 1:n_regions
        region_name = region_names{i};
        bounds = Constants.regions.(region_name);
        region_meas = data(bounds(1):bounds(2), :);
        
        T.(region_name) = (mean(region_meas, 1)).';
    end
end