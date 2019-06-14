%% Load Cata's data
i410_1_cata = csvread('/Users/sean/Desktop/test data cata analysis/i410_1.csv', 1, 1);
i410_1_cata = ssquare(i410_1_cata);

i410_2_cata = csvread('/Users/sean/Desktop/test data cata analysis/i410_2.csv', 1, 1);
i410_2_cata = ssquare(i410_2_cata);

abs_err_old = (abs(i410_1_cata - i410_2_cata));
%% Load my data
m_sep = readtable("sean/analysis/01-21-19 error analysis/data/movement_separated.csv");
rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
allImgs = loadTiffImageStack(rawImgFilePath);
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
im410_1_sean = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2_sean = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL"));

nAnimal1 = size(I.ImgFrame(I.LambdaGroup == "410/410" & I.Strain == "Animal 1" & I.SetFrame == 1), 1);
imTL = repmat(imTL(:,:,1), [1 1 nAnimal1]);

% Rotate Images according to first frame
seg410_1 = segmentPharynx(im410_1_sean, 0, 2000);

% im410_1 = rotatePharynx(im410_1, seg410_1);
% im410_2 = rotatePharynx(im410_2, seg410_1);
% imTL = rotatePharynx(imTL, seg410_1);

% Crop
% im410_1 = im410_1(52:79, 54:121, :); 
% im410_2 = im410_2(52:79, 54:121, :);
% imTL = imTL(52:79, 54:121, :);

% [i1, i2, matchingVecs, midlines1, midlines2, scaled_bounds1, scaled_bounds2, fdObjs, dx, dy, unreg_i1, unreg_i2] = pipelineTwoMidlinesTwoMasksRegistration(imTL, im410_1, im410_2);

[i410_1_2mid, i410_2_2mid, seg1, seg2, midlines1, midlines2, lrBounds1, lrBounds2, i1_raw, i2_raw] = pipelineTwoMidlinesTwoMasks(imTL, im410_1_sean, im410_2_sean);
[i410_1_sr, i410_2_sr, midlines1, midlines2, unreg_i1, unreg_i2, fdObjs] = pipelineSmoothRoughRegister(imTL, im410_1_sean, im410_2_sean);

% abs_err_2mid = (abs(i410_1_2mid - i410_2_2mid));%./((i410_1_2mid + i410_1_2mid)/2)) * 100;

% abs_err_sr = (abs(i410_1_sr - i410_2_sr));%./((i410_1_sr + i410_2_sr)/2)) * 100;

%% Configuring Plots

p_space_x_axis_label_size = 20;
intensity_y_axis_label_size = 20;

intensity_min = 0;
intensity_max = 130000;




%%
plotMultiplePharynxData({abs_err_old, (i410_1_cata + i410_2_cata)/2}, {'Abs. Error (Old Pipeline)', 'Avg. Intensity (Old Pipeline)'}, [0 13000]);

ylim_ = 13000;
ax = gca;
xlim = 100;

hXLabel = xlabel('P-A Axis of Pharynx');
hYLabel = ylabel('Intensity');
hLegend = findobj(gcf, 'Type', 'Legend');
yrule = ax.YAxis;
xrule = ax.XAxis;

yrule.FontSize = 20;
xrule.FontSize = 20;
set([hXLabel, hYLabel], 'FontSize', 22);
set(hLegend, 'FontSize', 22);
set([xrule, yrule, hLegend], 'FontWeight', 'bold');


set(gca, 'FontName', 'Helvetica');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
'XMinorTick', 'off', 'YMinorTick', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', round(linspace(0, ylim_, 5), 2), ...
'XTick', 0:25:xlim, 'LineWidth', 1);

%% Percent Error Plot

per_err_old = (abs(i410_1_cata - i410_2_cata)./((i410_1_cata + i410_2_cata)/2));

ylim_ = .20;
plotMultiplePharynxData({per_err_old}, {'Percent Error, Old Pipeline'}, [0 ylim_]);
ax = gca;
xlim = 100;

hXLabel = xlabel('P-A Axis of Pharynx');
hYLabel = ylabel('Percent Error');
hLegend = findobj(gcf, 'Type', 'Legend');
yrule = ax.YAxis;
xrule = ax.XAxis;

font_scale = 2;
yrule.FontSize = 20;
xrule.FontSize = 20;
set([hXLabel, hYLabel], 'FontSize', 22);
set(hLegend, 'FontSize', 22);
set([xrule, yrule, hLegend], 'FontWeight', 'bold');


set(gca, 'FontName', 'Helvetica');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
'XMinorTick', 'off', 'YMinorTick', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', round(linspace(0, ylim_, 5), 2), ...
'XTick', 0:25:100, 'LineWidth', 1);


%%
figure;
any_moving_idx = sum(table2array(m_sep) >= 2, 2) > 0;
post_0 = per_err_old(:,~any_moving_idx) * 100;
post_1 = per_err_old(:,any_moving_idx) * 100;
% post_2 = err_old(:,m_sep.PosteriorBulb == 2);
% post_3 = err_old(:,m_sep.PosteriorBulb == 3);

ylim_ = 20;
cmap = cbrewer('qual', 'Set1', 3);
plotMultiplePharynxData({post_0, post_1}, {"No movement", "Any movement"}, [0 ylim_], gca, cmap);

ax = gca;
xlim = 100;

hXLabel = xlabel('P-A Axis of Pharynx');
hYLabel = ylabel('Percent Error');
hLegend = findobj(gcf, 'Type', 'Legend');
yrule = ax.YAxis;
xrule = ax.XAxis;

font_scale = 2;
yrule.FontSize = 20;
xrule.FontSize = 20;
set([hXLabel, hYLabel], 'FontSize', 22);
set(hLegend, 'FontSize', 22);
set([xrule, yrule, hLegend], 'FontWeight', 'bold');


set(gca, 'FontName', 'Helvetica');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
'XMinorTick', 'off', 'YMinorTick', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', round(linspace(0, ylim_, 5), 2), ...
'XTick', 0:25:100, 'LineWidth', 1);

%%
abs_err_2mid = (abs(i410_1_2mid - i410_2_2mid));
abs_err_sr = (abs(i410_1_sr - i410_2_sr));

per_err_2mid = (abs(i410_1_2mid - i410_2_2mid)./((i410_1_2mid + i410_2_2mid)/2));
per_err_sr = (abs(i410_1_2mid - i410_2_2mid)./((i410_1_2mid + i410_2_2mid)/2));

errors_old = regionMeansLong(per_err_old, Constants.regions, 'PercentError', 'Old');
errors_2mid = regionMeansLong(per_err_2mid,  Constants.regions, 'PercentError', '2 Midlines');
errors_sr = regionMeansLong(per_err_sr,  Constants.regions, 'PercentError', 'SmoothRough');

big_table = vertcat(errors_old, errors_2mid, errors_sr);
writetable(big_table, '/Users/sean/Desktop/test data cata analysis/errors_natural.csv');

ylim_ = 20;
cmap = cbrewer('qual', 'Set1', 3);
plotMultiplePharynxData({100*per_err_old, 100*per_err_2mid}, {"Old Pipeline", "New Pipeline"}, [0 ylim_], gca, cmap);
ax = gca;
xlim = 100;

hXLabel = xlabel('P-A Axis of Pharynx');
hYLabel = ylabel('Percent Error');
hLegend = findobj(gcf, 'Type', 'Legend');
yrule = ax.YAxis;
xrule = ax.XAxis;

font_scale = 2;
yrule.FontSize = 20;
xrule.FontSize = 20;
set([hXLabel, hYLabel], 'FontSize', 22);
set(hLegend, 'FontSize', 22);
set([xrule, yrule, hLegend], 'FontWeight', 'bold');


set(gca, 'FontName', 'Helvetica');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
'XMinorTick', 'off', 'YMinorTick', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', round(linspace(0, ylim_, 5), 2), ...
'XTick', 0:25:100, 'LineWidth', 1);

%% Add registration
[reg_i1, reg_i2, ~, ~, fdObjs] = channelRegisterSmoothRough(i410_1_cata, i410_2_cata, 100);
[reg_i1_smooth, reg_i2_smooth, ~, ~, fdObjs_smooth] = ChannelRegister(i410_1_cata, i410_2_cata, 100);
%% Strategies by Movement

cmap = cbrewer('qual', 'Paired', 6);

i410_1_cata_reg = eval_fd(linspace(3,97,100), reg_i1);
i410_2_cata_reg = eval_fd(linspace(3,97,100), reg_i2);

abs_err_cata_reg = abs(i410_1_cata_reg - i410_2_cata_reg);
per_err_cata_reg = abs_err_cata_reg ./ ((i410_1_cata_reg + i410_2_cata_reg)/2);

err_no_movement_old  = 100 * per_err_old(:, ~any_moving_idx);
err_any_movement_old = 100 * per_err_old(:,  any_moving_idx);

err_no_movement_2mid  = 100 * per_err_2mid(:, ~any_moving_idx);
err_any_movement_2mid = 100 * per_err_2mid(:,  any_moving_idx);

err_no_movement_old_reg  = 100 * per_err_cata_reg(:, ~any_moving_idx);
err_any_movement_old_reg = 100 * per_err_cata_reg(:,  any_moving_idx);

ylim_ = 20;
plotMultiplePharynxData({
    err_no_movement_old,  err_any_movement_old,...
    err_no_movement_old_reg, err_any_movement_old_reg,...
    err_no_movement_2mid, err_any_movement_2mid},...
    {"Old (Stationary)", "Old (Moving)",...
    "Old with Registration (Stationary)", "Old with registration (Moving)",...
    "2 Midlines (Stationary)", "2 Midlines (Moving)"}, [0 ylim_], gca, cmap);


errors_old_stationary_tbl = regionMeansLong(err_no_movement_old, Constants.regions, 'PercentError', 'Old (Stationary)');
errors_old_moving_tbl     = regionMeansLong(err_any_movement_old, Constants.regions, 'PercentError', 'Old (Moving)');

errors_2mid_stationary_tbl = regionMeansLong(err_no_movement_2mid,  Constants.regions, 'PercentError', '2 Midlines (Stationary)');
errors_2mid_moving_tbl     = regionMeansLong(err_any_movement_2mid, Constants.regions, 'PercentError', '2 Midlines (Moving)');

errors_old_reg_stationary_tbl = regionMeansLong(err_no_movement_old_reg,  Constants.regions, 'PercentError', 'Old+Reg (Stationary)');
errors_old_reg_moving_tbl     = regionMeansLong(err_any_movement_old_reg, Constants.regions, 'PercentError', 'Old+Reg (Moving)');

big_table = vertcat(errors_old_stationary_tbl, errors_old_moving_tbl,...
    errors_2mid_stationary_tbl, errors_2mid_moving_tbl,...
    errors_old_reg_stationary_tbl, errors_old_reg_moving_tbl);

writetable(big_table, '/Users/sean/Desktop/test data cata analysis/errors_natural_by_movement.csv');

%% Smooth-rough demo
figure;
ylim_ = 14000;
hold on;
scatter(1:100, i410_1_cata(:,1), 'x');
plot(reg_i1(1)); hold on;
plot(reg_i1_smooth(1));
legend("Discrete Observations", "l=1", "l=100");
ylim([0 ylim_]);


cmap = cbrewer('qual', 'Set1', 3);
ax = gca;
xlim = 100;

hXLabel = xlabel('P-A Axis of Pharynx');
hYLabel = ylabel('Intensity');
hLegend = findobj(gcf, 'Type', 'Legend');
yrule = ax.YAxis;
xrule = ax.XAxis;

font_scale = 2;
yrule.FontSize = 20;
xrule.FontSize = 20;
set([hXLabel, hYLabel], 'FontSize', 22);
set(hLegend, 'FontSize', 22);
set([xrule, yrule, hLegend], 'FontWeight', 'bold');


set(gca, 'FontName', 'Helvetica');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
'XMinorTick', 'off', 'YMinorTick', 'off', ...
'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', round(linspace(0, ylim_, 5), 2), ...
'XTick', 0:25:100, 'LineWidth', 1);