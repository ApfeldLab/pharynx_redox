%% Load Data
p_410_m_410 = csvread('../experiments/single_v_dual_poly/poly_410_measure_410_sean.csv', 1, 1);
p_410_m_410(p_410_m_410<2000) = 0;

sq = flipud(ssquare(p_410_m_410));

subset = 1:10;
data = sq(:,subset);

%% Smooth

L = 0.0891; % Determined experimentally via findOptimalLambda.m
wormFd = makeWormFd_SJ(data, 'lambda', L);
% Unregistered
figure;
plot(wormFd);
title('Smoothed Functions (No Registration)');
drawnow;

%% Landmark Registration
landmarks = findLandmarks_SJ(wormFd);
% line([landmarks(idx,1) landmarks(idx,1)], get(gca, 'ylim'), 'Color', 'red', 'LineStyle','--');
% line([landmarks(idx,2) landmarks(idx,2)], get(gca, 'ylim'), 'Color', 'red', 'LineStyle','--');

% tic
[lregFd, lwarpFd] = landmarkreg(wormFd, landmarks);
% toc

% Landmarks
figure;
plot(lregFd);
title('Landmark-Registered Functions');
drawnow;

%% Continuous Registration

[cregFd, cwarpFd] = register_fd(mean(lregFd), lregFd, fdPar(getbasis(lregFd), int2Lfd(2), 200000));

% Continuous
figure;
plot(cregFd);
title('Continuous-Registered Functions (from Landmark Reg.)');
drawnow;
%% Plots

figure;
ax = subplot(1,1,1);
ylim(ax, [0 max(data(:))]);

for i=1:size(data,1)
    
    % Unregistered
    plot(ax, eval_fd(linspace(1,100,1000), mean(wormFd)), 'b--', 'LineWidth', 4); hold on;
    plot(ax, eval_fd(linspace(1,100,1000), wormFd(i)), 'b', 'LineWidth', 2); hold on;
    
    % Landmark Registered
    plot(ax, eval_fd(linspace(1,100,1000), mean(lregFd)), 'r--', 'LineWidth', 4); hold on;
    plot(ax, eval_fd(linspace(1,100,1000), lregFd(i)), 'r', 'LineWidth', 2); hold on;
    
    % Continuously Registered
    plot(ax, eval_fd(linspace(1,100,1000), mean(cregFd)), 'g--', 'LineWidth', 4); hold on;
    plot(ax, eval_fd(linspace(1,100,1000), cregFd(i)), 'g', 'LineWidth', 2); hold off;
    set(gca,'color',[0 0 0]);
    
    % Continuously Registered
    drawnow;
    pause;
end