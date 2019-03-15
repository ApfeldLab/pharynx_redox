%% FLAGS
LOAD_DATA = 0;
GENERATE_TABLES = 1;
WRITE_TABLES = 1;
REDO_CATA = 0;
REDO_1MID = 0;
REDO_2MID = 0;
REDO_REG  = 1;
PLOT = 0;

%% Load data
if LOAD_DATA
    [imTL,im410_1, im410_2, movement] = loadErrorData();
end
%% Run all pipelines
if REDO_CATA
    [i1_cata, i2_cata] = pipelineCata(imTL, im410_1, im410_2);
end

if REDO_1MID
    [i1_1mid, i2_1mid] = pipelineOneMidlineTwoMasks(imTL, im410_1, im410_2);
end

if REDO_2MID
    [i1_2mid, i2_2mid] = pipelineTwoMidlinesTwoMasks(imTL, im410_1, im410_2);
end

if REDO_REG
    [i1, i2, matchingVecs, mids1, mids2, scaled_bounds1, scaled_bounds2, fdObjs, dx, dy] = pipelineTwoMidlinesTwoMasksRegistration(imTL, im410_1, im410_2);
end

%% Tables
if GENERATE_TABLES
abs_err_cata = regionMeansLong(abs(i1_cata - i2_cata), Constants.regions, 'AbsoluteError', 'Cata');
abs_err_1mid = regionMeansLong(abs(i1_1mid - i2_1mid), Constants.regions, 'AbsoluteError', '1 Midline');
abs_err_2mid = regionMeansLong(abs(i1_2mid - i2_2mid), Constants.regions, 'AbsoluteError', '2 Midlines');
abs_err_reg  = regionMeansLong(abs(i1_reg - i2_reg),   Constants.regions, 'AbsoluteError', 'Registered');

abs_err_table = vertcat(abs_err_cata, abs_err_1mid, abs_err_2mid, abs_err_reg);

int_1_cata = regionMeansLong(i1_cata, Constants.regions, 'I1', 'Cata');
int_1_1mid = regionMeansLong(i1_1mid, Constants.regions, 'I1', '1 Midline');
int_1_2mid = regionMeansLong(i1_2mid, Constants.regions, 'I1', '2 Midlines');
int_1_reg  = regionMeansLong(i1_reg,  Constants.regions, 'I1', 'Registered');

int_1_table = vertcat(int_1_cata, int_1_1mid, int_1_2mid, int_1_reg);

int_2_cata = regionMeansLong(i2_cata, Constants.regions, 'I2', 'Cata');
int_2_1mid = regionMeansLong(i2_1mid, Constants.regions, 'I2', '1 Midline');
int_2_2mid = regionMeansLong(i2_2mid, Constants.regions, 'I2', '2 Midlines');
int_2_reg  = regionMeansLong(i2_reg,  Constants.regions, 'I2', 'Registered');

int_2_table = vertcat(int_2_cata, int_2_1mid, int_2_2mid, int_2_reg);

rel_err_cata = regionMeansLong(abs(i1_cata - i2_cata) ./ ((i1_cata + i2_cata) / 2), Constants.regions, 'RelativeError', 'Cata');
rel_err_1mid = regionMeansLong(abs(i1_1mid - i2_1mid) ./ ((i1_1mid + i2_1mid) / 2), Constants.regions, 'RelativeError', '1 Midline');
rel_err_2mid = regionMeansLong(abs(i1_2mid - i2_2mid) ./ ((i1_2mid + i2_2mid) / 2), Constants.regions, 'RelativeError', '2 Midlines');
rel_err_reg  = regionMeansLong(abs(i1_reg  - i2_reg)  ./ ((i1_reg  + i2_reg)  / 2), Constants.regions, 'RelativeError', 'Registered');

rel_err_table = vertcat(rel_err_cata, rel_err_1mid, rel_err_2mid, rel_err_reg);

big_table = join(abs_err_table, join(int_1_table, join(int_2_table, rel_err_table)));

% dist_means = regionMeansLong(d, Constants.regions, 'Distance', 'NA');
dx_means = regionMeansLong(dx, Constants.regions, 'dX', 'NA');
dy_means = regionMeansLong(dy, Constants.regions, 'dY', 'NA');

% big_table.Distance = repmat(dist_means.Distance, [4 1]);
big_table.dx = repmat(dx_means.dX, [4 1]);
big_table.dy = repmat(dy_means.dY, [4 1]);
end

if WRITE_TABLES
    writetable(big_table, '~/Desktop/normal_errors_table.csv');
end
%% Plot Intensity
if PLOT
data_ = cat(3,...
    i1_cata, i1_1mid, i1_2mid, i1_reg);

labels_ = {'Cata', '1 Midline, Two Masks (Unregistered)', '2 Midlines, Two Masks (Unregistered)', '2 Midlines, Two Masks (Registered)'};
ylim_ = [0 20000];
plotMultiplePharynxData(data_, labels_, ylim_);
addRegionBoundsToPlot(gca, Constants.regions);
end
%% Absolute Error
% data_ = cat(3,...
%     abs(i1_cata - i2_cata),...
%     abs(i1_1mid - i2_1mid),...
%     abs(i1_2mid - i2_2mid),...
%     abs(i1_reg - i2_reg));
% 
% labels_ = {'Cata', '1 Midline, Two Masks (Unregistered)', '2 Midlines, Two Masks (Unregistered)', '2 Midlines, Two Masks (Registered)'};
% ylim_ = [0 1000];
% plotMultiplePharynxData(data_, labels_, ylim_);
% addRegionBoundsToPlot(gca, Constants.regions);
