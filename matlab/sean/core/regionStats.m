function stats = regionStats(measurements, regions)
    fn = fieldnames(Constants.regions);
    stats = struct();
    for i=1:numel(fn)-2
        b = Constants.regions.(fn{i});
        
    end
end