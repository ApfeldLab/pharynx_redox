function HIDPlot(ratioImage, mask, im410_1, im410_2, ...
    i410_f1_m1, i410_f2_m1, i410_f2_m2,...
    i410_f1_cata, i410_f2_cata,...
    reg410_1, reg410_2, midline1, midline2,...
    scaledBounds1, scaledBounds2, warpFD, mvmt, i)
%HIDPLOT Summary of this function goes here
%   Detailed explanation goes here

cmap = cbrewer('qual', 'Dark2', 5);
set(groot, 'defaultAxesColorOrder', cmap);

% Ratio Image with Midlines
subplot(3, 6, [1 2]);
plotMultipleMidlines(ratioImage, im410_1, mask, midline1, midline2, scaledBounds1, scaledBounds2, warpFD);
mvmt = table2array(mvmt);
title(sprintf('[%d] | PB: %d | AB: %d | ST: %d | T: %d', i, mvmt(1), mvmt(2), mvmt(3), mvmt(4)));

% Intensities
subplot(3, 6, [3 6]);
cmap_paired = cbrewer('qual', 'Paired', 8);
set(groot, 'defaultAxesColorOrder', cmap_paired);
xs = 1:size(i410_f1_m1,1);
plot(xs, i410_f1_m1, xs, i410_f2_m2, xs, i410_f1_cata, xs, i410_f2_cata);
addRegionBoundsToPlot(gca, Constants.regions);
legend('F1 M1', 'F2 M2', 'F1 Cata', 'F2 Cata');
xlim([1 100]);
ylabel("I_{410}");


% Registered Intensities
set(groot, 'defaultAxesColorOrder', cmap);
subplot(3, 6, [9 12]);
xs = linspace(1,100,300);
hold on;
plot(xs, eval_fd(xs, reg410_1), xs, eval_fd(xs, reg410_2));
ylabel("I_{410}");
yyaxis right;
plot(xs, eval_fd(xs, warpFD));
ylim([1 100]);
xlim([1 100]);
addRegionBoundsToPlot(gca, Constants.regions);
hold off;
legend('r410_1', 'r410_2', 'warp', 'y=x');



% % Ratio Profiles
subplot(3,6, [15, 18]);
xs = 1:100;
% xss = linspace(1,100,300);

reg_ratios = eval_fd(xs, reg410_1 ./ reg410_2);

plot(xs, i410_f1_m1 ./ i410_f2_m1,...
    xs, i410_f1_m1 ./ i410_f2_m2,...
    xs, reg_ratios,...
    xs, i410_f1_cata ./ i410_f2_cata); hold on;
ylim([.8 1.2]);
hline(1);
addRegionBoundsToPlot(gca, Constants.regions);
legend('I1M1 / I2M1',  'I1M1 / I2M2', 'R1M1 / R2M2', 'I1C / I2C');
ylabel("I1 / I2");

%
subplot(3, 6, [13 14]);
single_mid_errors = regionMeans(abs(1 - (i410_f1_m1 ./ i410_f2_m1)));
double_mid_errors = regionMeans(abs(1 - (i410_f1_m1 ./ i410_f2_m2)));
reg_errors = regionMeans(abs(1 - reg_ratios));
cata_errors = regionMeans(abs(1 - (i410_f1_cata ./ i410_f2_cata)));

tbl = vertcat(single_mid_errors, double_mid_errors, reg_errors, cata_errors);
bar(table2array(tbl).');
set(gca,'xticklabel', fieldnames(Constants.regions));
legend('I1M1 / I2M1',  'I1M1 / I2M2', 'R1M1 / R2M2', 'I1C / I2C');
xtickangle(45);
ylabel("Avg. Error in Ratio");
end