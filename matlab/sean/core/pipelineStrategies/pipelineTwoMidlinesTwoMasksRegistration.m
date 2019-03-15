function [i1, i2, matchingVecs, midlines1, midlines2, scaled_bounds1, scaled_bounds2, fdObjs, dx, dy, unreg_i1, unreg_i2] = pipelineTwoMidlinesTwoMasksRegistration(imTL, imFL1, imFL2)
%PIPELINEONEMIDLINE Summary of this function goes here
%   Detailed explanation goes here
% Subtract Medians
imFL1 = subtractMedians(imFL1);
imFL2 = subtractMedians(imFL2);

% Segment
seg1 = segmentPharynx(imFL1, 0);
seg2 = segmentPharynx(imFL2, 0);

% Draw Midlines
disp('Calculating Midlines... Please hold');
midlines1 = calculateMidlines(imTL, seg1, imFL1, 0);
midlines2 = calculateMidlines(imTL, seg2, imFL2, 0);
disp('Done calculating midlines');

% Measure Under Midlines
N_DATA_POINTS = 1000;

lrBounds1 = double(getLeftRightBounds(seg1));
lrBounds2 = double(getLeftRightBounds(seg2));

[i1, ~, scaled_bounds1] = measureAndTrim(imFL1, midlines1, lrBounds1, N_DATA_POINTS);
[i2, ~, scaled_bounds2] = measureAndTrim(imFL2, midlines2, lrBounds2, N_DATA_POINTS);

unreg_i1 = i1;
unreg_i2 = i2;

[reg_i1, reg_i2, ~, ~, fdObjs] = ...
    ChannelRegister(i1, i2, 100);

i1 = eval_fd(1:100, reg_i1);
i2 = eval_fd(1:100, reg_i2);

matchingVecs = getMidlineMatchingVectors(midlines1, midlines2, scaled_bounds1, scaled_bounds2, fdObjs,100);
[dx, dy] = getDyDxMidlines(midlines1, midlines2, scaled_bounds1, scaled_bounds2, fdObjs,100);
end