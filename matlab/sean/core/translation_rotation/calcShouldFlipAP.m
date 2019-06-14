function [shouldFlip] = calcShouldFlipAP(intensityData)
%CALCSHOULDFLIPAP Summary of this function goes here
%   Detailed explanation goes here

profileLength = size(intensityData, 1);
nAnimals = size(intensityData, 2);
shouldFlip = false(nAnimals, 1);

for i=1:nAnimals
    minPeakDistance = 0.3 * size(intensityData, 1);
    [pks,locs_] = findpeaks(intensityData(:,i), 'MinPeakDistance', minPeakDistance);
    [~, pk_idx] = maxk(pks, 2);
    locs_ = locs_(pk_idx);
    if min(locs_) < profileLength - max(locs_)
        shouldFlip(i) = 1;
    end
end