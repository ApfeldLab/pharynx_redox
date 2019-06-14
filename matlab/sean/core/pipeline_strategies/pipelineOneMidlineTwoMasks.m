function [i1, i2] = pipelineOneMidlineTwoMasks(imTL, imFL1, imFL2)
%PIPELINEONEMIDLINE Summary of this function goes here
%   Detailed explanation goes here
% Subtract Medians
imFL1 = subtractMedians(imFL1);
imFL2 = subtractMedians(imFL2);

% Segment
seg1 = segmentPharynx(imFL1, 0);
seg2 = segmentPharynx(imFL2, 0);

% Draw Midlines
midlines1 = calculateMidlines(imTL, seg1, imFL1, 0);

% Measure Under Midlines
N_DATA_POINTS = 1000;

lrBounds1 = double(getLeftRightBounds(seg1));
lrBounds2 = double(getLeftRightBounds(seg2));

i1 = measureAndTrim(imFL1, midlines1, lrBounds1, N_DATA_POINTS);
i2 = measureAndTrim(imFL2, midlines1, lrBounds2, N_DATA_POINTS);
end

