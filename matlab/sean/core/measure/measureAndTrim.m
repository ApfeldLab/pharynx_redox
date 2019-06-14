function [i, i_raw] = measureAndTrim(imStack, midlines, trimThreshold, profileLength)
    i_raw = measureIntensityAlongMidlines(imStack, midlines, profileLength, 'BILINEAR');
    i_trimmed = trimProfile(i_raw, trimThreshold);
    i = ssquare(i_trimmed);
end