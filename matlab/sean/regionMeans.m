function T = regionMeans(d)
%REGIONMEANS Create a table, with headers for each region, and each row
%corresponding to the mean of the values of a single animal inside that
%region.
    region_names = fieldnames(Constants.regions);
    T = table;
    for i = 1:numel(region_names)
        region_name = region_names{i};
        region = Constants.regions.(region_name);

        T.(region_name) = (mean(d(region(1):region(2),:), 1)).';
    end
end