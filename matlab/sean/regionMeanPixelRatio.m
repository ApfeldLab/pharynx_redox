function T = regionMeanPixelRatio(intensityData)
%REGIONINTEGRATEDRATIO Create a table for the "Region-Mean Pixel Ratio",
%with a header for each region, each row being the n-th animal.
% 
%   To calculate the Mean Pixel Ratio, we divide each pixel in the 410nm
%   channel by the corresponding pixel in the 470nm channel. We then
%   calculate the mean of the resultant values.
% 
%   The Region-Mean Pixel Ratio is calculated along specific ranges. Those
%   regions are defined in the Constants class.
% 
%   ARGUMENTS:
%   	intensityData: a structure containing 2 fields: m410 and m470,
%   	representing measurements in the 410nm and 470nm wavelength
%   	channels. These should have identical dimensions (they should be
%       "squared" before processing by this function.

    fn = fieldnames(Constants.regions);
    T = table;
    for i = 1:numel(fn)
        region_name = fn{i};
        region = Constants.regions.(region_name);
        T.(strcat(region_name, '_mpbp')) = (mean(intensityData.m410(region(1):region(2),:) ./ intensityData.m470(region(1):region(2),:),1)).';
    end
end