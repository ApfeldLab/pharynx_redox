function [reg410_FD, reg470_FD, warp470, resampled_intensity] = ChannelRegister(...
    meas_410,meas_470, varargin)
%UNTITLED3 Register the 470 to the 410
%   Measurements should be the same dimensions (use square)

if (~isequal(size(meas_410), size(meas_470)))
    error('Measurement arrays must be the same dimensions');
end

% Register
WARP_NBASIS = 6;
WARP_ORDER = 4;
WARP_LAMBDA = 5000;

n_worms = size(meas_410, 2);

warpBasis = create_bspline_basis([1 100], WARP_NBASIS, WARP_ORDER);
fdParObj = fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA);

% TODO: Preallocate fdObjs for memory/speed

textprogressbar('Registering 470 to 410: ');
for i=1:n_worms
    textprogressbar(100*i/n_worms);
    origFd = makeWormFd_SJ(horzcat(meas_410(:,i),meas_470(:,i)), 'lambda', 10^0.0891);
    [regFd, warpFd] = register_fd(origFd(1), origFd, fdParObj, 0, 2, 1e-4, 100, 0);
    
    fdObjs(i).origFD = origFd;
    fdObjs(i).regFD = regFd;
    fdObjs(i).warpFD = warpFd;
end
textprogressbar('Done.');

% Resample
xs = linspace(1,100,100);
int_410 = zeros(size(xs, 2), n_worms);
int_470 = zeros(size(xs, 2), n_worms);
warp470 = zeros(size(xs, 2), n_worms);

for i=1:n_worms
    int_410(:,i) = eval_fd(xs, fdObjs(i).regFD(1));
    int_470(:,i) = eval_fd(xs, fdObjs(i).regFD(2));
    temp = eval_fd(xs, fdObjs(i).warpFD); % temp variable, this eval spits out two cols, only care about the 2nd
    warp470(:,i) = temp(:,2);
end

reg410_FD = makeWormFd_SJ(int_410, 'lambda', 10^0.0891);
reg470_FD = makeWormFd_SJ(int_470, 'lambda', 10^0.0891);

resampled_intensity.m410 = int_410;
resampled_intensity.m470 = int_470;

end