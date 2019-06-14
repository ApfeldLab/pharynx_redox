function [i1, i2, midlines1, midlines2, unreg_i1, unreg_i2, fdObjs] = pipelineSmoothRoughRegister(imTL, imFL1, imFL2)
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

[i1, ~, ~] = measureAndTrim(imFL1, midlines1, lrBounds1, N_DATA_POINTS);
[i2, ~, ~] = measureAndTrim(imFL2, midlines2, lrBounds2, N_DATA_POINTS);

unreg_i1 = i1;
unreg_i2 = i2;

smoothLambda = 2;
roughLambda = .0001;
warpLambda = 2;

[smoothFD1, smoothFD2, roughFD1, roughFD2, regRoughFD2, regSmoothFD2, warpFD] = smoothRoughRegister(i1, i2, smoothLambda, roughLambda, warpLambda);

fdObjs.smoothFD1 = smoothFD1;
fdObjs.smoothFD2 = smoothFD2;
fdObjs.roughFD1 = roughFD1;
fdObjs.roughFD2 = roughFD2;
fdObjs.regRoughFD2 = regRoughFD2;
fdObjs.regSmoothFD2 = regSmoothFD2;
fdObjs.warpFD = warpFD;

i1 = eval_fd(1:100, roughFD1);
i2 = eval_fd(1:100, regRoughFD2);

end