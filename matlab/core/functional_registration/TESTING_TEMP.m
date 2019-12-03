%% LOAD DATA

profile_data = ncread('/Users/sean/code/pharynx_redox/data/paired_ratio/2017_02_22-HD233_SAY47/analyses/2019-10-11_unreg/2017_02_22-HD233_SAY47-untrimmed_profile_data.nc', '__xarray_dataarray_variable__');

i410_0 = squeeze(profile_data(:, 1, 2, :));
i470_0 = squeeze(profile_data(:, 1, 1, :));
i410_1 = squeeze(profile_data(:, 2, 2, :));
i470_1 = squeeze(profile_data(:, 2, 1, :));

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

rgh_410_0 = makeWormFd_SJ(i410_0, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_breaks', ROUGH_N_BASIS);
% rgh_470_0 = makeWormFd_SJ(i470_0, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_basis', ROUGH_N_BASIS);
sm_410_0 = makeWormFd_SJ(i410_0, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BASIS);
% sm_470_0 = makeWormFd_SJ(i470_0, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_basis', SMOOTH_N_BASIS);
% 
% rgh_410_1 = makeWormFd_SJ(i410_1, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_basis', ROUGH_N_BASIS);
% rgh_470_1 = makeWormFd_SJ(i470_1, 'lambda', ROUGH_LAMBDA, 'n_order', ROUGH_ORDER, 'n_basis', ROUGH_N_BASIS);
% sm_410_1 = makeWormFd_SJ(i410_1, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_basis', SMOOTH_N_BASIS);
% sm_470_1 = makeWormFd_SJ(i470_1, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_basis', SMOOTH_N_BASIS);

figure;
to_register_idx = 1;
scatter(linspace(1, 100, 200), i410_0(:,to_register_idx)); hold on;
xs = linspace(1, 100, 1000);

plot(xs, eval_fd(xs, rgh_410_0(to_register_idx)));
plot(xs, eval_fd(xs, sm_410_0(to_register_idx)));

% figure;
% plot(xs, eval_fd(xs, deriv(rgh_410_0(to_register_idx), N_DERIV))); hold on;
% plot(xs, eval_fd(xs, deriv(rgh_470_0(to_register_idx), N_DERIV)));
% 
% title('rough');
% 
% figure;
% plot(xs, eval_fd(xs, deriv(sm_410_0(to_register_idx), N_DERIV)), '-r'); hold on;
% plot(xs, eval_fd(xs, deriv(sm_470_0(to_register_idx), N_DERIV)), '-.r');
% plot(xs, eval_fd(xs, deriv(sm_410_1(to_register_idx), N_DERIV)), '-b');
% plot(xs, eval_fd(xs, deriv(sm_470_1(to_register_idx), N_DERIV)), '-.b');
% title('smooth');
% 
% figure;
% plot(xs, eval_fd(xs, deriv(sm_410_0(to_register_idx), 0)), '-r'); hold on;
% plot(xs, eval_fd(xs, deriv(sm_470_0(to_register_idx), 0)), '-.r');
% plot(xs, eval_fd(xs, deriv(sm_410_1(to_register_idx), 0)), '-b');
% plot(xs, eval_fd(xs, deriv(sm_470_1(to_register_idx), 0)), '-.b');
% title('smooth');
% 
% 
% figure
% title('R');
% plot(xs, eval_fd(xs, sm_410_0(to_register_idx) ./ sm_470_0(to_register_idx)), '-r'); hold on;
% plot(xs, eval_fd(xs, sm_410_1(to_register_idx) ./ sm_470_1(to_register_idx)), '-b'); hold off;

%% REG ALL TO 410
to_register_idx = 11;
[preg_rgh_410_0, preg_rgh_410_1, pwarpfd_1, pwfd_1, pd410_1, pd470_1] = channel_register(...
    i410_0(:,to_register_idx), i410_1(:,to_register_idx), ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV);

[preg_rgh_410_0, preg_rgh_410_1, pwarpfd_1, pwfd_1, pd410_1, pd470_1] = channel_register(...
    i410_0(:,to_register_idx), i410_1(:,to_register_idx), ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV);


%% pair-register
% to_register_idx = 1:5;
to_register_idx = 11;
[preg_rgh_410_0, preg_rgh_410_1, pwarpfd_1, pwfd_1, pd410_1, pd470_1] = channel_register(...
    i410_0(:,to_register_idx), i410_1(:,to_register_idx), ...
    10, 4, 1,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    0);
%%
preg_rgh_470_0 = rgh_470_0;
preg_rgh_470_1 = synch(linspace(0, 1, 1000), rgh_470_1(to_register_idx), pwfd_1);

% i = 5 ;
figure;

% 410
plot(xs, eval_fd(xs, preg_rgh_410_0(1)), '-r'); hold on;
plot(xs, eval_fd(xs, rgh_410_1(to_register_idx)), '-b');
plot(xs, eval_fd(xs, preg_rgh_410_1(1)), '-.b');


% 470
figure;
plot(xs, eval_fd(xs, preg_rgh_470_0(to_register_idx)), '-r'); hold on;
plot(xs, eval_fd(xs, rgh_470_1(to_register_idx)), '-b');
plot(xs, eval_fd(xs, preg_rgh_470_1(1)), '-.b');

figure;
plot(xs, eval_fd(xs, rgh_410_0(to_register_idx) ./ rgh_470_0(to_register_idx)), '-r'); hold on;
plot(xs, eval_fd(xs, rgh_410_1(to_register_idx) ./ rgh_470_1(to_register_idx)), '-b');

plot(xs, eval_fd(xs, preg_rgh_410_0(1) ./ preg_rgh_470_0(to_register_idx)), '-.r');
plot(xs, eval_fd(xs, preg_rgh_410_1(1) ./ preg_rgh_470_1(1)), '-.b');

%% channel-register

[reg_rgh_410_0, reg_rgh_470_0, warpfd_0, wfd_0, d410_0, d470_0] = channel_register(...
    i410_0(:,to_register_idx), i470_0(:,to_register_idx), ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV);

[reg_rgh_410_1, reg_rgh_470_1, warpfd_1, wfd_1, d410_1, d470_1] = channel_register(...
    i410_1(:,to_register_idx), i470_1(:,to_register_idx), ...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA,...
    SMOOTH_LAMBDA, SMOOTH_N_BASIS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BASIS, ROUGH_ORDER,...
    N_DERIV);

%%
% xs = linspace(0, 1, 1000);
% scatter(linspace(0, 1, 200), i410_0(:,i)); hold on;
% plot(xs, eval_fd(xs, rgh_410_0(i)));
% 
% %%
% figure;
% plot(linspace(0,1,100),linspace(0,1,100)); hold on;
% plot(warpfd_0);
% 
% %%
% xs = linspace(0, 1, 1000);
% 
% i = 18;
% i = 10;
% figure;
% plot(xs, eval_fd(xs, d410_0(i))); hold on;
% plot(xs, eval_fd(xs, deriv(sm_470_0(i), N_DERIV)));
% plot(xs, eval_fd(xs, d470_0(i)));
% legend('d410', 'd470', 'r(470)');
% % 
% % %%
% xs = linspace(0, 1, 1000);
% figure;
% 
% scatter(linspace(0, 1, 200), i470_0(:, 1)); hold on;
% plot(xs, eval_fd(xs, reg_rgh_470_0(1))); 
% plot(xs, eval_fd(xs, rgh_470_0(1)));
% plot(xs, (eval_fd(xs, rgh_410_0(1))/1.48));
% legend('470', 'reg 470', 'raw 470', '410');
% 
% 
% %%
xs = linspace(0, 1, 1000);
figure;

to_register_idx = 1;
plot(xs, eval_fd(xs, reg_rgh_410_0(to_register_idx) ./ reg_rgh_470_0(to_register_idx)), '-.b'); hold on;
plot(xs, eval_fd(xs, reg_rgh_410_1(to_register_idx) ./ reg_rgh_470_1(to_register_idx)), '-.r');
plot(xs, eval_fd(xs, rgh_410_0(to_register_idx) ./ rgh_470_0(to_register_idx)),'-b');
plot(xs, eval_fd(xs, rgh_410_1(to_register_idx) ./ rgh_470_1(to_register_idx)),'-r');
legend('reg1', 'reg2', 'raw1', 'raw2');

%%
% xs = linspace(0, 1, 1000);
% figure;
% unreg_error = abs(eval_fd(xs, (rgh_410_0(to_register_idx) ./ rgh_470_0(to_register_idx)) - (rgh_410_1(to_register_idx) ./ rgh_470_1(to_register_idx))));
% reg_error = abs(eval_fd(xs, (rgh_410_0(to_register_idx) ./ rgh_470_0(to_register_idx)) - (reg_rgh_410_1(to_register_idx) ./ reg_rgh_470_1(to_register_idx))));
% plot(mean(unreg_error, 2)); hold on;
% plot(mean(reg_error, 2));
% legend('unreg', 'reg');

%%
% i = 11;
% plot(xs, eval_fd(xs, deriv(sm_410_1(i), 1))); hold on;
% plot(xs, eval_fd(xs, deriv(sm_470_1(i), 1)));


%%

xs = linspace(0, 1, 1000);
figure;

to_register_idx = 1;
% % unreg
% plot(xs, eval_fd(xs, (rgh_410_1(i) - mean(i410_1, 'all')) ./ std(i410_1, 0, 'all')), '-r');
% hold on;
% plot(xs, eval_fd(xs, (rgh_470_1(i) - mean(i470_1, 'all')) ./ std(i470_1, 0, 'all')), '-b');
% 
% % reg
% plot(xs, eval_fd(xs, (reg_rgh_410_1(i) - mean(i410_1, 'all')) ./ std(i410_1, 0, 'all')), '-.r');
% plot(xs, eval_fd(xs, (reg_rgh_470_1(i) - mean(i470_1, 'all')) ./ std(i470_1, 0, 'all')), '-.b'); 

% unreg
plot(xs, eval_fd(xs, rgh_410_1(11)), '-r'); hold on;
plot(xs, eval_fd(xs, rgh_470_1(11)), '-b');

% reg
plot(xs, eval_fd(xs, reg_rgh_410_1(to_register_idx)), '-.r');
plot(xs, eval_fd(xs, reg_rgh_470_1(to_register_idx)), '-.b');

% legend('410', '470', 'r410', 'r470');

hold off;

%%
to_register_idx = 1;
plot(rgh_410_0(to_register_idx) ./ rgh_470_0(to_register_idx)); hold on;
plot(rgh_410_1(to_register_idx) ./ rgh_470_1(to_register_idx)); 
legend('r0', 'r1');
hold off;

figure;
plot(reg_rgh_410_0(to_register_idx) ./ reg_rgh_470_0(to_register_idx)); hold on;
plot(reg_rgh_410_1(to_register_idx) ./ reg_rgh_470_1(to_register_idx)); 
legend('reg r0', 'reg r1');
hold off;

%%
% xs = linspace(0, 1, 1000);
% figure;
% unreg_r0 = rgh_410_0 ./ rgh_470_0;
% unreg_r1 = rgh_410_1 ./ rgh_470_1;
% 
% reg_r0 = reg_rgh_410_0 ./ reg_rgh_470_0;
% reg_r1 = reg_rgh_410_1 ./ reg_rgh_470_1;
% 
% plot(mean(unreg_r0)); hold on;
% plot(mean(unreg_r1));
% plot(mean(reg_r0));
% plot(mean(reg_r1));
% legend('ur0', 'ur1', 'r0', 'r1');
% ylim([1.4, 1.51]);
% xlim([10, 90]);