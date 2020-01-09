function [reg410_discrete, reg470_discrete, warp470_discrete] = ...
    channel_register(...
    meas_410, meas_470, ...
    resample_resolution, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA, ...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER, ...
    N_DERIV)

if (~isequal(size(meas_410), size(meas_470)))
    error('Measurement arrays must be the same dimensions');
end

n_worms = size(meas_410, 1);

% TODO make the basis range more general
warpBasis = create_bspline_basis([1 100], WARP_N_BASIS, WARP_ORDER);
warp_fdParObj = fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA);

disp('Registering 470 to 410: ');

[sm410, ~] = makeWormFd_SJ(meas_410, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
[sm470, ~] = makeWormFd_SJ(meas_470, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);

[rgh410, ~] = makeWormFd_SJ(meas_410, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
[rgh470, rgh470Par] = makeWormFd_SJ(meas_470, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);

[~, warpFd, wfd] = register_fd(deriv(sm410, N_DERIV), deriv(sm470, N_DERIV), warp_fdParObj);

putfd(rgh470Par, rgh470);

reg_rgh_470 = synch2(linspace(1, 100, 1000), rgh470Par, wfd);

% Resample
xs = linspace(1, 100, resample_resolution);
reg410_discrete = zeros(size(xs, 2), n_worms);
reg470_discrete = zeros(size(xs, 2), n_worms);
warp470_discrete = zeros(size(xs, 2), n_worms);

for i = 1:n_worms
    reg410_discrete(:, i) = eval_fd(xs, rgh410(i));
    reg470_discrete(:, i) = eval_fd(xs, reg_rgh_470(i));
    warp470_discrete(:, i) = eval_fd(xs, warpFd(i));
end

end
