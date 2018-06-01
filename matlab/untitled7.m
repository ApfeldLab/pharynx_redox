%% Mean over Regions
region_names = fieldnames(Constants.regions);
T = table;
d = data.E;
for i = 1:numel(region_names)
    region_name = region_names{i};
    region = Constants.regions.(region_name);
    
    T.(strcat('mean_',region_name)) = mean(d(region(1):region(2),:), 2);
end