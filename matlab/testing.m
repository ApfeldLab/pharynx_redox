
WARP_N_BASIS = 30;
WARP_ORDER = 4;
WARP_LAMBDA = 5e3;
SMOOTH_LAMBDA = 100;
SMOOTH_N_BREAKS = 100;
SMOOTH_ORDER = 4;
ROUGH_LAMBDA = 0.01;
ROUGH_N_BREAKS = 300;
ROUGH_ORDER = 4;
N_DERIV = 1;


i410 = squeeze(data(1:5, 1, 1, :));
i470 = squeeze(data(1:5, 1, 1, :));

[reg410_discrete, reg470_discrete, warp470_discrete] = channel_register(...
    i410, i470, 200, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA, ...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER, ...
    N_DERIV);