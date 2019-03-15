function intensityProfiles = measureIntensityCata(fluorescenceStack, maskStack, midlines, profileLength)
%MEASUREINTENSITYCATA Summary of this function goes here
%   Detailed explanation goes here
    masked = fluorescenceStack .* maskStack;
%     bounds = getLeftRightBounds(maskStack);
    bounds = repmat([1 size(fluorescenceStack, 2)], size(fluorescenceStack, 3), 1);
    intensityProfiles = measureIntensityAlongMidlines(masked, midlines, bounds, profileLength, 'BILINEAR');
    intensityProfiles(intensityProfiles<=2000) = 0;
    intensityProfiles = ssquare(intensityProfiles);
end