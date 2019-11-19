function [reg_rgh_410, reg_rgh_470, warpfd, wfd, d410, reg_d470] = ...
    channel_register_new(...
    meas_410, meas_470, ...
    warp_to, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER,...
    N_DERIV)

if (~isequal(size(meas_410), size(meas_470)))
    error('Measurement arrays must be the same dimensions');
end

warpBasis = create_bspline_basis([1 100], WARP_N_BASIS, WARP_ORDER);
fdParObj = fdPar(warpBasis, int2Lfd(2), WARP_LAMBDA);

disp('Registering 470 to 410: ');

sm410 = makeWormFd_SJ(meas_410, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
sm470 = makeWormFd_SJ(meas_470, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);

rgh410 = makeWormFd_SJ(meas_410, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
rgh470 = makeWormFd_SJ(meas_470, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);

d_target = deriv(warp_to, N_DERIV);

d410 = deriv(sm410, N_DERIV);
d470 = deriv(sm470, N_DERIV);

periodic = 0;
crit = 2;
conv = 1e-8;
iterlim = 50;
dbglev = 1;

target_warpBasis = create_bspline_basis([1 100], 4, 4);
target_fdParObj = fdPar(target_warpBasis, int2Lfd(2), 10);

disp('warping to target');
[t410, ~, wfd_target_410] = register_fd(d_target, d410, target_fdParObj, periodic, crit, conv, iterlim, dbglev);

disp('channel warping');
[reg_d470, warpfd, wfd] = register_fd(t410, d470, fdParObj, periodic, crit, conv, iterlim, dbglev);

reg_rgh_470 = synch(linspace(1,100,1000), rgh470, wfd);
reg_rgh_410 = synch(linspace(1,100,1000), rgh410, wfd_target_410);

end