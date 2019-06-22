function [i, i_raw] = measureAndTrim(imStack, midlines, trimThreshold, profileLength, measurementInterpMethod, resizeInterpMethod)
    i_raw = measureIntensityAlongMidlines(imStack, midlines, profileLength, measurementInterpMethod);
    i_trimmed = trimProfile(i_raw, trimThreshold);
    i = ssquare(i_trimmed, profileLength, resizeInterpMethod);
end