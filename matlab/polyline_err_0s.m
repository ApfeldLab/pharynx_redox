coords_410 = csvread('../experiments/single_v_dual_poly/sean_coords_410.csv', 1, 1);
coords_470 = csvread('../experiments/single_v_dual_poly/sean_coords_470.csv', 1, 1);
labels = csvread('../data/raw_data_1/labels.csv');

%%

c_410_x = coords_410(:, 1:2:end);
c_410_y = coords_410(:, 2:2:end);

c_470_x = coords_410(:, 1:2:end);
c_470_y = coords_470(:, 2:2:end);

%%
errs = c_410_y - c_470_y;
plot(abs(mean(errs(:,labels==0), 2))); hold on;

plot(abs(mean(errs(:,labels==1), 2))); hold on;

plot(abs(mean(errs(:,labels==2), 2))); hold on;

plot(abs(mean(errs(:,labels==3), 2)));

legend('0', '1', '2', '3');

%%
plot(errs(:,labels==0), 'Color', 'b'); hold on;
plot(errs(:,labels==1), 'Color', 'm', 'LineWidth', 2); hold on;
plot(errs(:,labels==2), 'Color', 'c', 'LineWidth', 2); hold on;
plot(errs(:,labels==3), 'Color', 'g', 'LineWidth', 2);