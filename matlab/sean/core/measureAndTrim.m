function [i, unscaled_bounds, scaled_bounds, i_raw] = measureAndTrim(imStack, midlines, lrBounds, profileLength)
    i_raw = measureIntensityAlongMidlines(imStack, midlines, lrBounds, profileLength, 'BILINEAR');
    [i_trimmed, unscaled_bounds, scaled_bounds] = trimProfile(i_raw, lrBounds);
    i = ssquare(i_trimmed);
end