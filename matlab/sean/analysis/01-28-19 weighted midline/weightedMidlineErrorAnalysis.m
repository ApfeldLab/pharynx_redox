% This script runs a pipeline for the following imaging approach:
%
% A SINGLE transmitted light image is captured for each animal. A 
% series of alternating 410nm and 470nm images are taken of this animal.
% These are treated as pairs. Next, a series of alternating 410nm and 410nm
% images are taken of the same animal. Again, these are treated as pairs.
%
% The experimenter must provide the indices at which a new animal starts.
% This index metadata is contained in a csv file with the following format
% 
%   start:end,ID,lambda1/lambda2/lambda3 (etc)
%   
% So, for example, the file might look like this:
%       1:1,     Animal 1, TL
%       2:127,   Animal 1, 470/410a
%       128:241, Animal 1, 410/410
%       242:242, Animal 2, TL
%       243:364, Animal 2, 470/410
%       365:488, Animal 2, 410/410
%
%
% The start:end column MUST contain a start AND end, even if it is a
% single frame. In that case, use the same frame for start and end (see
% example 1, rows 1 and 4).

analysis_dir = '/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-28-19 weighted midline';

% Load images
rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-28-19 weighted midline/data/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);

allImgs = loadTiffImageStack(rawImgFilePath);

% Load Movement
m_sep = readtable("sean/analysis/01-21-19 error analysis/data/movement_separated.csv");

% Split Images
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
[im410_animal1_frame1_raw, im410_animal1_frame2_raw, imTL_animal1] = splitImages(allImgs, "Animal 1", I);
[im410_animal2_frame1_raw, im410_animal2_frame2_raw, imTL_animal2] = splitImages(allImgs, "Animal 2", I);

nAnimals = size(im410_animal1_frame1_raw, 3);

% Subtract Medians
im410_animal1_frame1 = subtractMedians(im410_animal1_frame1_raw);
im410_animal1_frame2 = subtractMedians(im410_animal1_frame2_raw);
im410_animal2_frame1 = subtractMedians(im410_animal2_frame1_raw);
im410_animal2_frame2 = subtractMedians(im410_animal2_frame2_raw);

% Segment
seg410_animal1_frame1 = segmentPharynx(im410_animal1_frame1, 0);
seg410_animal1_frame2 = segmentPharynx(im410_animal1_frame2, 0);
seg410_animal2_frame1 = segmentPharynx(im410_animal2_frame1, 0);
seg410_animal2_frame2 = segmentPharynx(im410_animal2_frame2, 0);
%%
% Draw Midlines
segLRBounds_animal1_frame1 = double(getLeftRightBounds(seg410_animal1_frame1));
segLRBounds_animal1_frame2 = double(getLeftRightBounds(seg410_animal1_frame2));
segLRBounds_animal2_frame1 = double(getLeftRightBounds(seg410_animal2_frame1));
segLRBounds_animal2_frame2 = double(getLeftRightBounds(seg410_animal2_frame2));

midlines410_animal1_frame1 = calculateMidlines(imTL_animal1, seg410_animal1_frame1, im410_animal1_frame1, 0);
midlines410_animal1_frame2 = calculateMidlines(imTL_animal1, seg410_animal1_frame2, im410_animal1_frame2, 0);
midlines410_animal2_frame1 = calculateMidlines(imTL_animal2, seg410_animal2_frame1, im410_animal2_frame1, 0);
midlines410_animal2_frame2 = calculateMidlines(imTL_animal2, seg410_animal2_frame2, im410_animal2_frame2, 0);

midlines410_animal1_frame1_weighted = calculateMidlines(imTL_animal1, seg410_animal1_frame1, im410_animal1_frame1, 1);
midlines410_animal1_frame2_weighted = calculateMidlines(imTL_animal1, seg410_animal1_frame2, im410_animal1_frame2, 1);
midlines410_animal2_frame1_weighted = calculateMidlines(imTL_animal2, seg410_animal2_frame1, im410_animal2_frame1, 1);
midlines410_animal2_frame2_weighted = calculateMidlines(imTL_animal2, seg410_animal2_frame2, im410_animal2_frame2, 1);

% Measure Under Midlines
N_DATA_POINTS = 1000;

[i410_animal_1_frame_1, ~, scaledLrBounds_animal_1_frame_1] = measureAndTrim(im410_animal1_frame1, midlines410_animal1_frame1, segLRBounds_animal1_frame1, N_DATA_POINTS);
[i410_animal_1_frame_2, ~, scaledLrBounds_animal_1_frame_2] = measureAndTrim(im410_animal1_frame2, midlines410_animal1_frame2, segLRBounds_animal1_frame2, N_DATA_POINTS);
[i410_animal_1_frame_2_midline_1, ~, scaledLrBounds_animal_1_frame_2_midline_1] = measureAndTrim(im410_animal1_frame2, midlines410_animal1_frame1, segLRBounds_animal1_frame1, N_DATA_POINTS);


[i410_animal_2_frame_1, ~, scaledLrBounds_animal_2_frame_1] = measureAndTrim(im410_animal2_frame1, midlines410_animal2_frame1, segLRBounds_animal2_frame1, N_DATA_POINTS);
[i410_animal_2_frame_2, ~, scaledLrBounds_animal_2_frame_2] = measureAndTrim(im410_animal2_frame2, midlines410_animal2_frame2, segLRBounds_animal2_frame2, N_DATA_POINTS);


[i410_animal_1_frame_1_weighted_midlines, ~, scaledLrBounds_animal_1_frame_1_weighted_midlines] = measureAndTrim(im410_animal1_frame1, midlines410_animal1_frame1_weighted, segLRBounds_animal1_frame1, N_DATA_POINTS);
[i410_animal_1_frame_2_weighted_midlines, ~, scaledLrBounds_animal_1_frame_2_weighted_midlines] = measureAndTrim(im410_animal1_frame2, midlines410_animal1_frame2_weighted, segLRBounds_animal1_frame2, N_DATA_POINTS);

i_animal_1_frame_1_cata = measureIntensityCata(im410_animal1_frame1, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
i_animal_1_frame_2_cata = measureIntensityCata(im410_animal1_frame2, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);

%%
% Register
[reg_i410_animal1_frame1, reg_i410_animal1_frame2, ~, ~, fdObjs_animal1] = ...
    ChannelRegister(i410_animal_1_frame_1, i410_animal_1_frame_2, 100);
[reg_i410_animal1_frame1_mid1, reg_i410_animal1_frame2_mid1, ~, ~, fdObjs_animal1_single_midline] = ...
    ChannelRegister(i410_animal_1_frame_1, i410_animal_1_frame_2_midline_1, 100);
[reg_i410_animal2_frame1, reg_i410_animal2_frame2, ~, ~, fdObjs_animal2] = ...
    ChannelRegister(i410_animal_2_frame_1, i410_animal_2_frame_2, 100);


[reg_i410_animal1_frame1_weighted_midlines, reg_i410_animal1_frame2_weighted_midlines, ~, ~, fdObjs_animal1_weighted_midlines] = ...
    ChannelRegister(i410_animal_1_frame_1_weighted_midlines, i410_animal_1_frame_2_weighted_midlines, 100);

%%
% Transform to Redox Potential
unreg_ratio_animal_1 = i410_animal_1_frame_1 ./ i410_animal_1_frame_2;
unreg_ratio_animal_1_single_midline = i410_animal_1_frame_1 ./ i410_animal_1_frame_2_midline_1;

unreg_ratio_animal_2 = i410_animal_2_frame_1 ./ i410_animal_2_frame_2;

reg_ratio_animal_1 = reg_i410_animal1_frame1 ./ reg_i410_animal1_frame2;
reg_ratio_animal_2 = reg_i410_animal2_frame1 ./ reg_i410_animal2_frame2;

%% Individual Error
unreg_error_animal_1 = i410_animal_1_frame_1 - i410_animal_1_frame_2;
unreg_error_animal_1_single_midline = i410_animal_1_frame_1 - i410_animal_1_frame_2_midline_1;
reg_error_animal_1 = eval_fd(1:100, reg_i410_animal1_frame1 - reg_i410_animal1_frame2);
reg_error_animal_1_single_midline = eval_fd(1:100, reg_i410_animal1_frame1_mid1 - reg_i410_animal1_frame2_mid1);
reg_error_animal_1_weighted_double_midline_registered = eval_fd(1:100, reg_i410_animal1_frame1_weighted_midlines - reg_i410_animal1_frame2_weighted_midlines);

registration_error_animal_1_frame_1 = eval_fd(1:100, reg_i410_animal1_frame1) - i410_animal_1_frame_1;


% Calculate Region Statistics
i410_animal1_frame1_regions = regionMeans(i410_animal_1_frame_1, Constants.regions, 'i410_1_');
i410_animal1_frame2_regions = regionMeans(i410_animal_1_frame_2, Constants.regions, 'i410_2_');

i410_animal1_frame1_regions_weighted = regionMeans(i410_animal_1_frame_1, Constants.regions, 'i410_1_weighted_');
i410_animal1_frame2_regions_weighted = regionMeans(i410_animal_1_frame_2, Constants.regions, 'i410_2_weighted_');

unreg_error_regions = regionMeans(unreg_error_animal_1, Constants.regions, 'unreg_error_');
reg_error_regions = regionMeans(reg_error_animal_1, Constants.regions, 'reg_error_');
unreg_error_single_midline_regions = regionMeans(unreg_error_animal_1_single_midline, Constants.regions, 'unreg_error_single_mid_');
reg_error_animal_1_single_midline_regions = regionMeans(reg_error_animal_1_single_midline, Constants.regions, 'reg_error_single_mid_');

reg_error_animal_1_weighted_double_midline_regions = regionMeans(reg_error_animal_1_weighted_double_midline_registered, Constants.regions, 'reg_error_double_weighted_');

% Distance
midline_dist = calcMidlineDistance(midlines410_animal1_frame1, midlines410_animal1_frame2, scaledLrBounds_animal_1_frame_1, scaledLrBounds_animal_1_frame_2, fdObjs_animal1);
midline_dist_regions = regionMeans(midline_dist, Constants.regions, 'midline_dist_');

errors_regions_table = horzcat(m_sep, i410_animal1_frame1_regions, unreg_error_regions, reg_error_regions, unreg_error_single_midline_regions, reg_error_animal_1_single_midline_regions, midline_dist_regions);

midline_strategies = {...
    'US', 'UD', ...
    'RS', 'RD', ...
    'WD', 'FN'
};

midline_strategies_errors = struct(...
    'US', regionMeans(unreg_error_animal_1_single_midline),...
    'UD', regionMeans(unreg_error_animal_1),...
    'RS', regionMeans(reg_error_animal_1_single_midline),...
    'RD', regionMeans(reg_error_animal_1), ...
    'WD', regionMeans(reg_error_animal_1_weighted_double_midline_registered),...
    'FN', regionMeans(registration_error_animal_1_frame_1));

midline_strategies_i410_1 = struct(...
    'US', regionMeans(i410_animal_1_frame_1),...
    'UD', regionMeans(i410_animal_1_frame_1),...
    'RS', regionMeans(eval_fd(1:100, reg_i410_animal1_frame1_mid1)),...
    'RD', regionMeans(eval_fd(1:100, reg_i410_animal1_frame1)),...
    'WD', regionMeans(eval_fd(1:100, reg_i410_animal1_frame1_weighted_midlines)),...
    'FN', regionMeans(eval_fd(1:100, reg_i410_animal1_frame1)));

midline_strategies_i410_2 = struct(...
    'US', regionMeans(i410_animal_1_frame_2_midline_1),...
    'UD', regionMeans(i410_animal_1_frame_2),...
    'RS', regionMeans(eval_fd(1:100, reg_i410_animal1_frame2_mid1)),...
    'RD', regionMeans(eval_fd(1:100, reg_i410_animal1_frame2)),...
    'WD', regionMeans(eval_fd(1:100, reg_i410_animal1_frame2_weighted_midlines)),...
    'FN', regionMeans(i410_animal_1_frame_1));

region_names = fieldnames(Constants.regions);
mvmt_region_names = m_sep.Properties.VariableNames;

nRows = size(i410_animal1_frame1_regions,1) * length(region_names) * length(mvmt_region_names) * length(midline_strategies);
%% Write Table
long_table = table('Size', [nRows 6], ...
    'VariableTypes', {'uint16', 'string', 'string', 'double', 'double', 'double'},...
    'VariableNames', {'pair_id', 'region', 'midline_strategy', 'error', 'i410_1', 'i410_2'});
count = 0;
for i = 1:size(i410_animal1_frame1_regions,1)
    for j = 1:length(region_names)
        for k = 1:length(midline_strategies)
            count = count + 1;
            error_  = midline_strategies_errors.(midline_strategies{k}).(region_names{j})(i);
            i410_1_ = midline_strategies_i410_1.(midline_strategies{k}).(region_names{j})(i);
            i410_2_ = midline_strategies_i410_2.(midline_strategies{k}).(region_names{j})(i);

            long_table.pair_id(count) = i;
            long_table.region(count) = region_names(j);
            long_table.midline_strategy(count) = midline_strategies(k);
            long_table.error(count) = error_;
            long_table.i410_1(count) = i410_1_;
            long_table.i410_2(count) = i410_2_;
        end
    end
end

writetable(errors_regions_table, fullfile(analysis_dir, 'analysis', 'errors_regions_table.csv'));
writetable(long_table, fullfile(analysis_dir, 'analysis', 'long_errors_regions_table.csv'));

%% Metadata table

return;
%% Generate Figures
figure;
ax = subplot(3,1,1);
plotPharynxDataByMovement(unreg_error_animal_1, m_sep.PosteriorBulb, ax);
title(ax, 'Unregistered Error in Animal 1');
ylim(ax, [0 500]);

ax = subplot(3,1,2);
plotPharynxDataByMovement(reg_error_animal_1, m_sep.PosteriorBulb, ax);
title(ax, 'Registered Error in Animal 1')
ylim(ax, [0 500]);

ax = subplot(3,1,3);
plotPharynxDataByMovement(reg_error_animal_1_single_midline, m_sep.PosteriorBulb, ax);
title(ax, 'Registered Error in Animal 1 with Single Midline');
ylim(ax, [0 500]);

%% 
plotPharynxDataSubsets(cat(3, ...
    i410_animal_1_frame_1,...
    i410_animal_1_frame_2, ...
    i410_animal_1_frame_2_midline_1, ...
    eval_fd(1:100, reg_i410_animal1_frame2)));
legend('I_1', 'I_2', 'I_2 (mid_1)', 'I_2_{reg}');

plotPharynxDataSubsets(cat(3,...
    eval_fd(1:100, reg_i410_animal1_frame1),...
    eval_fd(1:100, reg_i410_animal1_frame2)));
legend('Reg. I_1', 'Reg. I_2');

plotPharynxDataSubsets(cat(3, ...
    i410_animal_1_frame_1,...
    i410_animal_1_frame_2));
legend('I_1', 'I_2');


figure;
ax = subplot(2, 1, 1);
plotPharynxDataSubsets(cat(3, ...
    unreg_error_animal_1,...
    unreg_error_animal_1_single_midline), ax);
legend('E', 'E (single Midline)');

ax = subplot(2, 1, 2);
plotPharynxDataSubsets(cat(3, ...
    i410_animal_1_frame_1,...
    i410_animal_1_frame_2, ...
    i410_animal_1_frame_2_midline_1), ax);
legend('I_1', 'I_2 (mid_2)', 'I_2 (mid_1)');

%%
plotPharynxDataSubsets(cat(3, ...
    i410_animal_1_frame_1,...
    eval_fd(1:100,reg_i410_animal1_frame1)));
legend('I_1', 'f(I_1)');