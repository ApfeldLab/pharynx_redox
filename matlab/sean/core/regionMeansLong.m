function T = regionMeansLong(data, regions, label, strategy)
%REGIONMEANS Create a table, with headers for each region, and each row
%corresponding to the mean of the values of a single animal inside that
%region.
    if exist('regions', 'var')
        regions_ = regions;
    else
        regions_ = Constants.regions;
    end
    
    region_names = fieldnames(regions_);
    nAnimals = size(data,2);
    T = table('Size', [(nAnimals * length(region_names)) 4], 'VariableTypes', {'int16', 'double', 'string', 'string'}, 'VariableNames', {'Animal', label, 'Strategy', 'Region'});
    row = 1;
    for i = 1:nAnimals
        for j = 1:length(region_names)
            region_name = region_names{j};
            region_bounds = regions_.(region_name);

            val = mean(data(region_bounds(1):region_bounds(2), i), 1);

            T(row, 1) = {i};
            T(row, 2) = {val};
            T(row, 3) = {strategy};
            T(row, 4) = {region_name};
            row = row + 1;
        end
    end
end