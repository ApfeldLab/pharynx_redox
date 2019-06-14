function [i1, i2, seg1, seg2, midlines1, midlines2, lrBounds1, lrBounds2, i1_raw, i2_raw] = pipelineTwoMidlinesTwoMasks(imTL, imFL1, imFL2)
%PIPELINEONEMIDLINE Summary of this function goes here
%   Detailed explanation goes here
% Subtract Medians
imFL1 = subtractMedians(imFL1);
imFL2 = subtractMedians(imFL2);

% Segment
seg1 = segmentPharynx(imFL1, 0);
seg2 = segmentPharynx(imFL2, 0);

% Rotate
[rot_imFL1, rot_seg1] = rotatePharynx(imFL1, seg1);
[rot_imFL2, rot_seg2] = rotatePharynx(imFL2, seg2);

% Draw Midlines
disp('Calculating Midlines... Please hold');
% midlines1 = calculateMidlines(imTL, seg1, imFL1, 0);
% midlines2 = calculateMidlines(imTL, seg2, imFL2, 0);

midlines1 = calculateMidlinesNoTL(rot_seg1);
% midlines2 = calculateMidlinesNoTL(rot_seg2);
midlines2 = midlines1;
disp('Done calculating midlines');

% Measure Under Midlines
N_DATA_POINTS = 1000;

lrBounds1 = double(getLeftRightBounds(seg1));
lrBounds2 = double(getLeftRightBounds(seg2));

[i1, i1_raw] = measureAndTrim(rot_imFL1, midlines1, 2000, N_DATA_POINTS);
[i2, i2_raw] = measureAndTrim(rot_imFL2, midlines2, 2000, N_DATA_POINTS);
end