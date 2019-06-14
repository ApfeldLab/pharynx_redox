function a = j_resample(series, step, col)
%function a = j_resample(series, step, col)
%Resamples data as ceil(series/step)*step.
%ceil(A) rounds the elements of A to the nearest integers >= A.
%If col = 0 or > #cols in series, or is left out, then all cols
%are resampled; otherwise just series(:,col) is resampled.
%If step and col are left out then step = 1.

a = [];
if nargin < 2
    step=1;
end
if nargin < 3 | col == 0 | col>size(series,2)
    a = ceil(series/step)*step;
else
    a = series;
    a(:,col) = ceil(series(:,col)/step)*step;
end
