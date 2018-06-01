function T = regionIntegratedRatio(intensityData)
%REGIONINTEGRATEDRATIO Create a table for the "Region-Integrated Ratio",
%with a header for each region, each row being the n-th animal.
% 
%   The Integrated Ratio is defined by the scalar: SUM(410_xi) / SUM(470_xi).
%   Then the Region-Integrated Ratios are calculated along specific ranges
%   for xi. Those regions are defined in the Constants class.
% 
%   ARGUMENTS:
%   	intensityData: a structure containing 2 fields: m410 and m470,
%   	representing measurements in the 410nm and 470nm wavelength
%   	channels. These should have identical dimensions (they should be
%       "squared" before processing by this function.

    region_names = fieldnames(Constants.regions);
    T = table;
    for i = 1:numel(region_names)
        region_name = region_names{i};
        region = Constants.regions.(region_name);
        T.(strcat('R',region_name)) = (sum(intensityData.m410(region(1):region(2),:)) ./ sum(intensityData.m470(region(1):region(2),:))).';
    end
end

