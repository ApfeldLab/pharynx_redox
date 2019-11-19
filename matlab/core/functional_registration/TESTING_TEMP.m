%% LOAD DATA

profile_data = ncread('/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/analyses/2019-10-11_unreg/2017_02_22-HD233_SAY47-untrimmed_profile_data.nc', '__xarray_dataarray_variable__');

i410_0 = squeeze(profile_data(:, 1, 2, :));
i470_0 = squeeze(profile_data(:, 1, 1, :));
i410_1 = squeeze(profile_data(:, 2, 2, :));
i470_1 = squeeze(profile_data(:, 2, 1, :));

r_0 = i410_0 ./ i470_0;
r_1 = i410_1 ./ i470_1;

%% Registration testing


% figure;
% ax = gca;

% plot(ax, ...
%     xs, i410_0(:,1),...
%     xs, i470_0(:,1),...
%     xs, eval_fd(xs, rgh_410_0(1)),...
%     xs, eval_fd(xs, sm_410_0(1))...
% );

% plot(ax, xs, eval_fd(xs, sm_410_0(1), 1)); hold on;
% plot(ax, xs, eval_fd(xs, sm_470_0(1), 1));

WARP_N_BASIS = 10;
WARP_ORDER = 4;
WARP_LAMBDA = 5;

SMOOTH_LAMBDA = 2;
SMOOTH_N_BREAKS = 50;
SMOOTH_ORDER = 4;

ROUGH_LAMBDA = .1;
ROUGH_N_BREAKS = 200;
ROUGH_ORDER = 4;

N_DERIV = 1;

rgh_410_0 = makeWormFd_SJ(i410_0, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
rgh_470_0 = makeWormFd_SJ(i470_0, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
sm_410_0 = makeWormFd_SJ(i410_0, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
sm_470_0 = makeWormFd_SJ(i470_0, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);

rgh_410_1 = makeWormFd_SJ(i410_1, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
rgh_470_1 = makeWormFd_SJ(i470_1, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BREAKS);
sm_410_1 = makeWormFd_SJ(i410_1, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
sm_470_1 = makeWormFd_SJ(i470_1, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);

i = 1;

figure;
scatter(linspace(1, 100, 200), i410_0(:,i)); hold on;
plot(xs, eval_fd(xs, rgh_410_0(i)));
plot(xs, eval_fd(xs, sm_410_0(i)));

figure;
plot(xs, eval_fd(xs, deriv(rgh_410_0(i), N_DERIV))); hold on;
plot(xs, eval_fd(xs, deriv(rgh_470_0(i), N_DERIV)));
title('rough');

figure;
plot(xs, eval_fd(xs, deriv(sm_410_0(i), N_DERIV))); hold on;
plot(xs, eval_fd(xs, deriv(sm_470_0(i), N_DERIV)));
title('smooth');

%% register
to_register_idx = 1:123;

warp_to = mean(rgh_410_0);

[reg_rgh_410_0, reg_rgh_470_0, warpfd_0, wfd_0, d410_0, d470_0] = channel_register_new(...
    i410_0(:,to_register_idx), i470_0(:,to_register_idx), ...
    warp_to, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER,...
    N_DERIV);

[reg_rgh_410_1, reg_rgh_470_1, warpfd_1, wfd_1, d410_1, d470_1] = channel_register_new(...
    i410_1(:,to_register_idx), i470_1(:,to_register_idx), ...
    warp_to, ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER,...
    N_DERIV);

%%
xs = linspace(1, 100, 1000);
scatter(linspace(1, 100, 200), i410_0(:,i)); hold on;
plot(xs, eval_fd(xs, rgh_410_0(i)));

%%
figure;
plot(linspace(1,100,100),linspace(1,100,100)); hold on;
plot(warpfd_0);

%%
xs = linspace(1, 100, 1000);

% i = 18;
i = 2;
figure;
plot(xs, eval_fd(xs, d410_0(i))); hold on;
plot(xs, eval_fd(xs, deriv(sm_470_0(i), N_DERIV)));
plot(xs, eval_fd(xs, d470_0(i)));

legend('410', '470', 'r(470)');

%%
xs = linspace(1, 100, 1000);
figure;

scatter(linspace(1, 100, 200), i470_0(:, 1)); hold on;
plot(xs, eval_fd(xs, reg_rgh_470_0(1))); 
plot(xs, eval_fd(xs, rgh_470_0(1)));
legend('reg 470', 'raw 470');


%%
xs = linspace(1, 100, 1000);
figure;

i = 1;
plot(xs, eval_fd(xs, reg_rgh_410_0(i) ./ reg_rgh_470_0(i)), '-.b'); hold on;
plot(xs, eval_fd(xs, reg_rgh_410_1(i) ./ reg_rgh_470_1(i)), '-.r');
plot(xs, eval_fd(xs, rgh_410_0(i) ./ rgh_470_0(i)),'-b');
plot(xs, eval_fd(xs, rgh_410_1(i) ./ rgh_470_1(i)),'-r');
legend('reg1', 'reg2', 'raw1', 'raw2');

%%
xs = linspace(1, 100, 1000);
figure;
unreg_error = abs(eval_fd(xs, (rgh_410_0 ./ rgh_470_0) - (rgh_410_1 ./ rgh_470_1)));
reg_error = abs(eval_fd(xs, (reg_rgh_410_0 ./ reg_rgh_470_0) - (reg_rgh_410_1 ./ reg_rgh_470_1)));
plot(mean(unreg_error, 2));  hold on;
plot(mean(reg_error, 2));
legend('unreg', 'reg');