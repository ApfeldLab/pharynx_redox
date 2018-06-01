function plotRegionBoundaries(regions)
%PLOTREGIONBOUNDARIES Plots lines where the region boundaries are 
%   regions is a struct, with name->[start end]
    fn = fieldnames(regions);
    colors = ['b', 'k', 'r', 'g', 'c'];
    for k=1:numel(fn)
        if ~isequal(fn{k}, 'medial_axis')
            c = colors(mod(k,numel(colors))+1);
            label = fn{k};
            vline(regions.(label), c, label);
        end
    end
end