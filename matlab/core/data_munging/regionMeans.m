function T = regionMeans(data, regions, prefix)
%REGIONMEANS Create a table, with headers for each region, and each row
%corresponding to the mean of the values of a single animal inside that
%region.
    if exist('regions', 'var')
        regions_ = regions;
    else
        regions_ = Constants.regions;
    end
    
    region_names = fieldnames(regions_);
    T = table;
    for i = 1:numel(region_names)
        region_name = region_names{i};
        region = regions_.(region_name);

        T.(region_name) = (mean(data(region(1):region(2),:), 1)).';
    end
    
    if exist('prefix', 'var')
        newRegionNames = region_names;
        for i=1:length(region_names)
            newRegionNames{i} = strcat(prefix, newRegionNames{i});
        end
        T.Properties.VariableNames = newRegionNames;
    end
end