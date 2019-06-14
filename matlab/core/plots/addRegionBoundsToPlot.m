function addRegionBoundsToPlot(ax,regions)
%ADDREGIONBOUNDSTOPLOT Summary of this function goes here
%   Detailed explanation goes here
leg = findobj(gcf, 'Type', 'Legend');
leg.AutoUpdate = 'off';
fields = fieldnames(regions);
nRegions = size(fields, 1) - 2;
cmap = cbrewer('qual', 'Dark2', nRegions);
    for i=1:nRegions
        region = regions.(fields{i});
        xline(ax, region(1), '--', fields{i}, 'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'bottom', 'Color', cmap(i, :));
        xline(ax, region(2), '--', 'Color', cmap(i, :));
    end
leg.AutoUpdate = 'on';
end