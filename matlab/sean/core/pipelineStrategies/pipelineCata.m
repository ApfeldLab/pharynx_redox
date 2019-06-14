function [i1,i2] = pipelineCata(imTL, imFL1, imFL2)
%PIPELINECATA Summary of this function goes here
%   Assuming the the images are in TL/410/410

% Subtract Medians
imFL1 = subtractMedians(imFL1);
imFL2 = subtractMedians(imFL2);

% Segment
seg = segmentPharynx(imFL1, 0, 2000);

% Mask
imFL1 = imFL1 .* seg;
imFL2 = imFL2 .* seg;

% Rotate
[rot_FL1, ~] = rotatePharynx(imFL1, seg);
[rot_FL2, rot_seg] = rotatePharynx(imFL2, seg);

% Draw Midlines
% midlines = calculateMidlines(imTL, seg, imFL1, 0);
midlines = calculateMidlinesNoTL(rot_seg);

% Measure Under Midlines
N_DATA_POINTS = 1000;

i1 = measureIntensityCata(rot_FL1, rot_seg, midlines, N_DATA_POINTS);
i2 = measureIntensityCata(rot_FL2, rot_seg, midlines, N_DATA_POINTS);
end