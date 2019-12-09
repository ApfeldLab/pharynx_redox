function [reg410_discrete, reg470_discrete, warp470_discrete] = ...
    channel_register_standard(...
    meas_410, meas_470, ...
    resample_resolution,...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER,...
    N_DERIV)

if (~isequal(size(meas_410), size(meas_470)))
    error('Measurement arrays must be the same dimensions');
end

n_worms = size(meas_410, 2);

warpBasis = create_bspline_basis([1 100], WARP_N_BASIS, WARP_ORDER);
fdParObj = fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA);

disp('Registering 470 to 410: ');

fd410 = makeWormFd_SJ(meas_410, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
fd470 = makeWormFd_SJ(meas_470, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);

[r470, warpFd, wfd] = register_fd(deriv(fd410, N_DERIV), deriv(fd470, N_DERIV), fdParObj);

% r470 = synch(linspace(1, 100, 1000), fd470, wfd);

% Resample
xs = linspace(1,100,resample_resolution);
reg410_discrete = zeros(size(xs, 2), n_worms);
reg470_discrete = zeros(size(xs, 2), n_worms);
warp470_discrete = zeros(size(xs, 2), n_worms);

for i=1:n_worms
    reg410_discrete(:,i) = eval_fd(xs, fd410(i));
    reg470_discrete(:,i) = eval_fd(xs, r470(i));
    warp470_discrete(:,i) = eval_fd(xs, warpFd(i));
end

end