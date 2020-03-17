WARP_LAMBDA = 5;
WARP_N_BASIS = 10;
WARP_ORDER = 4;

SMOOTH_LAMBDA = 1;
SMOOTH_N_BREAKS = 30;
SMOOTH_ORDER = 4;

ROUGH_LAMBDA = 0.001;
ROUGH_N_BREAKS = 200;
ROUGH_ORDER = 4;

N_DERIV = 2;

% [sm410, ~] = makeWormFd_SJ(i410, 'lambda', SMOOTH_LAMBDA, 'n_order', SMOOTH_ORDER, 'n_breaks', SMOOTH_N_BREAKS);
% plot(sm410);

%%
i410 = squeeze(data(:, 1, 1, 1, :));
i470 = squeeze(data(:, 1, 1, 2, :));

%%
[reg410_discrete, reg470_discrete, warp_discrete] = pop_register(...
    i410, i470, ...
    1000,...
    WARP_N_BASIS, WARP_ORDER, WARP_LAMBDA, ...
    SMOOTH_LAMBDA, SMOOTH_N_BREAKS, SMOOTH_ORDER, ...
    ROUGH_LAMBDA, ROUGH_N_BREAKS, ROUGH_ORDER, ...
    N_DERIV);

%%
figure; plot(i410'); title('raw410') ; figure ; plot(i470'); title('raw470')
figure; plot(reg410_discrete); title('reg410') ; figure; plot(reg470_discrete); title('reg470');

%%
figure;
plot(linspace(0, 1, 2e2), mean(i410, 1)); hold on;
plot(linspace(0, 1, 1e3), mean(reg410_discrete, 2));
legend('raw', 'reg');

%%
figure;
plot(linspace(0, 1, 2e2), mean(i410./i470, 1)); hold on;
plot(linspace(0, 1, 1e3), mean(reg410_discrete./reg470_discrete, 2));
legend('raw', 'reg');

%%
figure;
plot((i410./i470)');
title('raw');
figure;
plot(reg410_discrete./reg470_discrete);
title('reg');
% legend('raw', 'reg');






