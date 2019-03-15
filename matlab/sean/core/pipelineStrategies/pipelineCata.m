function [i1,i2] = pipelineCata(imTL, imFL1, imFL2)
%PIPELINECATA Summary of this function goes here
%   Assuming the the images are in TL/410/410

% Subtract Medians
imFL1 = subtractMedians(imFL1);
imFL2 = subtractMedians(imFL2);

% Segment
seg = segmentPharynx(imFL1, 0, 2000);

% Draw Midlines
midlines = calculateMidlines(imTL, seg, imFL1, 0);

% Measure Under Midlines
N_DATA_POINTS = 1000;


i1 = measureIntensityCata(imFL1, seg, midlines, N_DATA_POINTS);
i2 = measureIntensityCata(imFL2, seg, midlines, N_DATA_POINTS);
end

