%% LOAD DATA

% profile_data = ncread('/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/analyses/2019-10-11_unreg/2017_02_22-HD233_SAY47-untrimmed_profile_data.nc', '__xarray_dataarray_variable__');

i410_0 = squeeze(prof_raw(:, 1, 2, :));
i470_0 = squeeze(prof_raw(:, 1, 1, :));
i410_1 = squeeze(prof_raw(:, 2, 2, :));
i470_1 = squeeze(prof_raw(:, 2, 1, :));

r_0 = i410_0 ./ i470_0;
r_1 = i410_1 ./ i470_1;

%% Registration testing

WARP_N_BASIS = 300;
WARP_ORDER = 4;
WARP_LAMBDA = 1e5;

SMOOTH_LAMBDA = 1e2;
SMOOTH_N_BASIS = 100;
SMOOTH_ORDER = 4;

ROUGH_LAMBDA = 10^0.05;
ROUGH_N_BASIS = 300;
ROUGH_ORDER = 4;

N_DERIV = 2;

to_register_idx = 1;

%% REG ALL TO 410
to_register_idx = 11;
[reg410_discrete, reg470_discrete, warp470_discrete] = channel_register(...
    i410_0(to_register_idx, :), i410_1(to_register_idx, :), ...
    200, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV);