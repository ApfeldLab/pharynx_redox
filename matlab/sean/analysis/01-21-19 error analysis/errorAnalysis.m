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
%       2:127,   Animal 1, 470/410
%       128:241, Animal 1, 410/410
%       242:242, Animal 2, TL
%       243:364, Animal 2, 470/410
%       365:488, Animal 2, 410/410
%
%
% The start:end column MUST contain a start AND end, even if it is a
% single frame. In that case, use the same frame for start and end (see
% example 1, rows 1 and 4).

% Load images
rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);

allImgs = loadTiffImageStack(rawImgFilePath);

% Load Movement
m_sep = readtable("sean/analysis/01-21-19 error analysis/data/movement_separated.csv");


% Split Images
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
[im410_animal1_frame1_raw, im410_animal1_frame2_raw, imTL_animal1] = splitImages(allImgs, "Animal 1", I);
[im410_animal2_frame1_raw, im410_animal2_frame2_raw, imTL_animal2] = splitImages(allImgs, "Animal 2", I);

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

midlines410_animal1_frame1 = calculateMidlines(imTL_animal1, seg410_animal1_frame1, im410_animal1_frame1);
midlines410_animal1_frame2 = calculateMidlines(imTL_animal1, seg410_animal1_frame2, im410_animal1_frame2);
midlines410_animal2_frame1 = calculateMidlines(imTL_animal2, seg410_animal2_frame1, im410_animal2_frame1);
midlines410_animal2_frame2 = calculateMidlines(imTL_animal2, seg410_animal2_frame2, im410_animal2_frame2);

% Measure Under Midlines
N_DATA_POINTS = 1000;

[i410_animal_1_frame_1, ~, scaledLrBounds_animal_1_frame_1] = measureAndTrim(im410_animal1_frame1, midlines410_animal1_frame1, segLRBounds_animal1_frame1, N_DATA_POINTS);
[i410_animal_1_frame_2, ~, scaledLrBounds_animal_1_frame_2] = measureAndTrim(im410_animal1_frame2, midlines410_animal1_frame2, segLRBounds_animal1_frame2, N_DATA_POINTS);
[i410_animal_1_frame_2_midline_1, ~, scaledLrBounds_animal_1_frame_2_midline_1] = measureAndTrim(im410_animal1_frame2, midlines410_animal1_frame1, segLRBounds_animal1_frame1, N_DATA_POINTS);


[i410_animal_2_frame_1, ~, scaledLrBounds_animal_2_frame_1] = measureAndTrim(im410_animal2_frame1, midlines410_animal2_frame1, segLRBounds_animal2_frame1, N_DATA_POINTS);
[i410_animal_2_frame_2, ~, scaledLrBounds_animal_2_frame_2] = measureAndTrim(im410_animal2_frame2, midlines410_animal2_frame2, segLRBounds_animal2_frame2, N_DATA_POINTS);

%%
% Register
[reg_i410_animal1_frame1, reg_i410_animal1_frame2, ~, ~, fdObjs_animal1] = ...
    ChannelRegister(i410_animal_1_frame_1, i410_animal_1_frame_2, 100);
[reg_i410_animal1_frame1_mid1, reg_i410_animal1_frame2_mid1, ~, ~, fdObjs_animal1_single_midline] = ...
    ChannelRegister(i410_animal_1_frame_1, i410_animal_1_frame_2_midline_1, 100);
[reg_i410_animal2_frame1, reg_i410_animal2_frame2, ~, ~, fdObjs_animal2] = ...
    ChannelRegister(i410_animal_2_frame_1, i410_animal_2_frame_2, 100);


%%
% Transform to Redox Potential
unreg_ratio_animal_1 = i410_animal_1_frame_1 ./ i410_animal_1_frame_2;
unreg_ratio_animal_1_single_midline = i410_animal_1_frame_1 ./ i410_animal_1_frame_2_midline_1;

unreg_ratio_animal_2 = i410_animal_2_frame_1 ./ i410_animal_2_frame_2;

reg_ratio_animal_1 = reg_i410_animal1_frame1 ./ reg_i410_animal1_frame2;
reg_ratio_animal_2 = reg_i410_animal2_frame1 ./ reg_i410_animal2_frame2;

%% Individual Error
unreg_error_animal_1 = abs(i410_animal_1_frame_1 - i410_animal_1_frame_2);
unreg_error_animal_1_single_midline = abs(i410_animal_1_frame_1 - i410_animal_1_frame_2_midline_1);
reg_error_animal_1 = abs(eval_fd(1:100, reg_i410_animal1_frame1 - reg_i410_animal1_frame2));
reg_error_animal_1_single_midline = abs(eval_fd(1:100, reg_i410_animal1_frame1_mid1 - reg_i410_animal1_frame2_mid1));


% Calculate Region Statistics
i410_animal1_frame1_regions = regionMeans(i410_animal_1_frame_1, Constants.regions, 'i410_1_');
i410_animal1_frame2_regions = regionMeans(i410_animal_1_frame_2, Constants.regions, 'i410_2_');
unreg_error_regions = regionMeans(unreg_error_animal_1, Constants.regions, 'unreg_error_');
reg_error_regions = regionMeans(reg_error_animal_1, Constants.regions, 'reg_error_');
unreg_error_single_midline_regions = regionMeans(unreg_error_animal_1_single_midline, Constants.regions, 'unreg_error_single_mid_');
reg_error_animal_1_single_midline_regions = regionMeans(reg_error_animal_1_single_midline, Constants.regions, 'reg_error_single_mid_');

% Distance
midline_dist = calcMidlineDistance(midlines410_animal1_frame1, midlines410_animal1_frame2, scaledLrBounds_animal_1_frame_1, scaledLrBounds_animal_1_frame_2, fdObjs_animal1);
midline_dist_regions = regionMeans(midline_dist, Constants.regions, 'midline_dist_');

%
analysis_dir = 'sean/analysis/01-21-19 error analysis/analysis';
% writetable(i410_animal1_frame1_regions, fullfile(analysis_dir, 'i410_regions.csv'));
% writetable(unreg_error_regions, fullfile(analysis_dir, 'unreg_error_regions.csv'));
% writetable(reg_error_regions, fullfile(analysis_dir, 'reg_error_regions.csv'));

errors_regions_table = horzcat(m_sep, i410_animal1_frame1_regions, unreg_error_regions, reg_error_regions, unreg_error_single_midline_regions, reg_error_animal_1_single_midline_regions, midline_dist_regions);

%%
% long_table = table('VariableNames', {'animal', 'region', 'midline_strategy', 'mvmt_region', 'mvmt', 'error', 'midline_dist'});

midline_strategies = {...
    'US', 'UD', ...
    'RS', 'RD'
};

midline_strategies_errors = struct(...
    'US', regionMeans(unreg_error_animal_1_single_midline, Constants.regions, ''),...
    'UD', regionMeans(unreg_error_animal_1, Constants.regions, ''),...
    'RS', regionMeans(reg_error_animal_1_single_midline, Constants.regions, ''),...
    'RD', regionMeans(reg_error_animal_1, Constants.regions, '')...
);

midline_strategies_i410_1 = struct(...
    'US', regionMeans(i410_animal_1_frame_1)
    'UD',...
    'RS',...
    'RD',...
);

region_names = fieldnames(Constants.regions);
mvmt_region_names = m_sep.Properties.VariableNames;


nRows = size(i410_animal1_frame1_regions,1) * length(region_names) * length(mvmt_region_names) * length(midline_strategies);

long_table = table('Size', [nRows 7], ...
    'VariableTypes', {'uint16', 'string', 'string', 'string', 'uint8', 'double', 'double'},...
    'VariableNames', {'pair_id', 'region', 'midline_strategy', 'movement_region', 'movement', 'error', 'i410_1'});
long_cell = cell(nRows, 6);
count = 0;
for i = 1:size(i410_animal1_frame1_regions,1)
    for j = 1:length(region_names)
        for k = 1:length(mvmt_region_names)
            for l = 1:length(midline_strategies)
                count = count + 1;
                
                mvmnt_ = m_sep.(mvmt_region_names{k})(i);
                error_ = midline_strategies_errors.(midline_strategies{l}).(region_names{j})(i);
                i410_1_ = i410_animal1_frame1_regions.(strcat('i410_1_', region_names{j}))(i);
                i410_2_ = i410_animal1_frame2_regions.(strcat('i410_2_', region_names{j}))(i);
                
                
                long_table.pair_id(count) = i;
                long_table.region(count) = region_names(j);
                long_table.midline_strategy(count) = midline_strategies(l);
                long_table.movement_region(count) = mvmt_region_names(k);
                long_table.movement(count) = mvmnt_;
                long_table.error(count) = error_;
                long_table.i410_1(count) = i410_1_;
                long_table.i410_2(count) = i410_2_;
            end
        end
    end
end


writetable(errors_regions_table, fullfile(analysis_dir, 'errors_regions_table.csv'));
writetable(long_table, fullfile(analysis_dir, 'long_errors_regions_table.csv'));

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

%%

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPERS

%% `Info` Plot

% Show image with colormap
imR_animal_1 = im410_animal1_frame1_raw ./ im410_animal1_frame2_raw;


cmap = cbrewer('div', 'RdBu', 1000);
% cmap_scaled = MC_nonlinearCmap(cmap, 1, [.5 1.5], 1, 1);
[muHat,sigmaHat,muCI,sigmaCI] = normfit(imR_animal_1(:));
low = muHat - 2.576 * sigmaHat;
hi = muHat + 2.576 * sigmaHat;
figure;
set(gcf, 'Position', get(0, 'Screensize'));

for i=1:size(imR_animal_1, 3)
    clf;
    
    % R(s)
    subplot(2,3,1:2);
    hold on;
    plot(unreg_ratio_animal_1(:,i));
    plot(unreg_ratio_animal_1_single_midline(:,i), 'Color', 'black');
    ylim([0 2]);
    legend('1st/1st', '1st/2nd');
    title('R_{410/410}(s) for Midline Strategies');
    
    % I(s)
    subplot(2,3,4:5);
    hold on;
    plot(i410_animal_1_frame_1(:,i));
    plot(i410_animal_1_frame_2_midline_1(:,i));
    plot(i410_animal_1_frame_2(:,i));
    ylim([0 15000]);
    xlim([0 100]);
    legend('1st', '2nd (mid 1)', '2nd (mid 2)');
    title('I_{410}(s)');
    
    % R(x,y)
    subplot(2,3,3);
    I = imR_animal_1(:,:,i);
    I_masked = I.*seg410_animal1_frame1(:,:,i);
    I_masked(I_masked == 0) = NaN;
    imAlpha=ones(size(I_masked));
    imAlpha(isnan(I_masked))=0;

    bbox = regionprops(seg410_animal1_frame1(:,:,1), 'BoundingBox');
    bbox = floor(bbox.BoundingBox);

    buffer = 15;

    left = bbox(1) - buffer;
    top = bbox(2) - buffer;
    right = bbox(1) + bbox(3) + buffer;
    bottom = bbox(2) + bbox(4) + buffer;
    
    hold on;
%     imagesc(I_masked(top:bottom,left:right), 'AlphaData', imAlpha(top:bottom,left:right));
    imagesc(I_masked);
    
    % Plot Midlines
    xs410 = segLRBounds_animal1_frame1(i, 1):segLRBounds_animal1_frame1(i,2);
    ys410 = feval(midlines410_animal1_frame1{i}, xs410);
    
    plot(xs410, ys410);
    
    
%     plot(midlines410_animal1_frame2{i});
    hold off;
    
    colormap(cmap); colorbar;
    caxis([low hi]);
    set(gca,'color',0*[1 1 1]);
    axis square;
    title('R_{410/410}(x,y)');
    
    % Overall Figure
    sgtitle(sprintf('Animal %d', i));
    
    % Save
%     export_fig('/Users/sean/Desktop/meeting_12_6_18/figures/ratio_images/info.pdf', '-opengl', '-append');
    waitforbuttonpress;
end

function plotPharynxDataByMovement(r, mvmt_mode, ax)
    x=1:100;
    ts = tinv([0.025  0.975],size(r, 1)-1);

    y_data = zeros(100, 3);
    SEM_data = zeros(100, 3);
    CI_data = zeros(100, 2, 3);
    e_data = zeros(100, 2, 3);
    
    for i=1:4
        idx = mvmt_mode == i - 1;
        y_data(:,i) = mean(r(:,idx), 2).';
        SEM_data(:,i) = std(r(:,idx), 0, 2)/sqrt(size(r, 1));
        CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
        e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    
    cmap_scaled = flipud(cbrewer('qual', 'Dark2', 4));
    if nargin > 2
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled, ax);
        xlim(ax, [1 100]);
    else
        figure;
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled);
        xlim([1 100]);
    end
    
    legend('0', '1', '2', '3', 'Location', 'northwest');
    xlabel('Pharynx Space');
end

function plotPharynxDataSubsets(data3d, ax)
%     Each group should be a new 2d plane of pharynx data, 
%       (Worms x Pharynx x Group)

    x=1:100;
    ts = tinv([0.025  0.975],size(data3d, 1)-1);

    y_data = zeros(100, 3);
    SEM_data = zeros(100, 3);
    CI_data = zeros(100, 2, 3);
    e_data = zeros(100, 2, 3);
    
    for i=1:size(data3d, 3)
        y_data(:,i) = mean(data3d(:,:,i), 2).';
        SEM_data(:,i) = std(data3d(:,:,i), 0, 2)/sqrt(size(data3d, 1));
        CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
        e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    
    cmap_scaled = flipud(cbrewer('qual', 'Dark2', max(3, size(data3d, 3))));
    
    if nargin > 1
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled, ax);
        xlim(ax, [1 100]);
    else
        figure;
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled);
        xlim([1 100]);
    end
    
    
end

function [im410_1, im410_2, imTL] = splitImages(allImages, group_id_str, indexer)
    % NOTE: The way these images are split are specific to this experiment!
    % Read the documentation on the loadIndexer method for more info on
    % splitting up images with the indexer.
    %
    % For this particular analysis, we are interested in the 410/410 data set.
    % I am also splitting up this data on a per-animal basis.
    im410_1 = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "410/410" & indexer.SetFrame == 1 & indexer.Strain == group_id_str));
    im410_2 = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "410/410" & indexer.SetFrame == 2 & indexer.Strain == group_id_str));

    imTL = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "TL" & indexer.Strain == group_id_str));
    imTL = repmat(imTL, [1 1 size(im410_1, 3)]);
end

function [i, unscaled_bounds, scaled_bounds] = measureAndTrim(imStack, midlines, lrBounds, profileLength)
    i_raw = measureIntensityAlongMidlines(imStack, midlines, lrBounds, profileLength, 'BILINEAR');
    [i_trimmed, unscaled_bounds, scaled_bounds] = trimProfile(i_raw, lrBounds);
    i = ssquare(i_trimmed);
end