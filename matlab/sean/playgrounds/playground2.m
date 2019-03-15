%% Watershed
bw = seg(:,:,2);
D = bwdist(~bw);
figure
imshow(D,[],'InitialMagnification','fit')
title('Distance transform of ~bw')
D = -D;
D(~bw) = Inf;
L = watershed(D, 8);
L(~bw) = 0;
rgb = label2rgb(L,'jet',[.5 .5 .5]);
figure
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of D')

%% Load data
rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
allImgs = loadTiffImageStack(rawImgFilePath);
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL"));
im410_1 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));

nAnimal1 = size(I.ImgFrame(I.LambdaGroup == "410/410" & I.Strain == "Animal 1" & I.SetFrame == 1), 1);
imTL = repmat(imTL(:,:,1), [1 1 nAnimal1]);

movement = readtable("~/Desktop/meeting_12_6_18/data/movement_separated.csv");
%% Make New Image Stacks from All Permutation Pairs
pairsAnimal1 = unorderedPairs(1:nAnimal1);

im410_1_dupes = im410_1(:,:,pairsAnimal1(:,1));
im410_2_dupes = im410_1(:,:,pairsAnimal1(:,2));
imTL_dupes = repmat(imTL(:,:,1), [1 1 size(im410_1_dupes, 3)]);

imR = im410_1_dupes ./ im410_2_dupes;

%% Cata
[i1_cata, i2_cata] = pipelineCata(imTL_dupes, im410_1_dupes, im410_2_dupes);

% 1 Midline, Two Masks, No Registration
[i1_1mid, i2_1mid] = pipelineOneMidlineTwoMasks(imTL_dupes, im410_1_dupes, im410_2_dupes);

% 2 Midlines, Two Masks, No Registration
[i1_2mid, i2_2mid] = pipelineTwoMidlinesTwoMasks(imTL_dupes, im410_1_dupes, im410_2_dupes);

% 2 Midlines, Two Masks, Registration
[i1_reg, i2_reg, d_dupes] = pipelineTwoMidlinesTwoMasksRegistration(imTL_dupes, im410_1_dupes, im410_2_dupes);

% Table Dupes
err_cata_dupes = regionMeansLong(i1_cata - i2_cata, Constants.regions, 'Error', 'Cata');
err_1mid_dupes = regionMeansLong(i1_1mid - i2_1mid, Constants.regions, 'Error', '1 Midlines');
err_2mid_dupes = regionMeansLong(i1_2mid - i2_2mid, Constants.regions, 'Error', '2 Midlines');
err_reg_dupes = regionMeansLong(i1_reg - i2_reg, Constants.regions, 'Error', 'Registered');

% dist_means_dupes = regionMeansLong(d_dupes, Constants.regions, 'Distance', 'NA');
%
big_table_dupes = vertcat(err_cata_dupes, err_1mid_dupes, err_2mid_dupes, err_reg_dupes);
% big_table_dupes.Distance = repmat(dist_means_dupes.Distance, [4 1]);
writetable(big_table_dupes, '~/Desktop/dupes_errors_table.csv');
%% Non-Duped
[i1_cata_normal, i2_cata_normal] = pipelineCata(imTL, im410_1, im410_2);
[i1_1mid_normal, i2_1mid_normal] = pipelineOneMidlineTwoMasks(imTL, im410_1, im410_2);
[i1_2mid_normal, i2_2mid_normal] = pipelineTwoMidlinesTwoMasks(imTL, im410_1, im410_2);
[i1_reg_normal, i2_reg_normal, d_normal] = pipelineTwoMidlinesTwoMasksRegistration(imTL, im410_1, im410_2);

%% Tables for normal
err_cata_normal = regionMeansLong(i1_cata_normal - i2_cata_normal, Constants.regions, 'Error', 'Cata');
err_1mid_normal = regionMeansLong(i1_1mid_normal - i2_1mid_normal, Constants.regions, 'Error', '1 Midlines');
err_2mid_normal = regionMeansLong(i1_2mid_normal - i2_2mid_normal, Constants.regions, 'Error', '2 Midlines');
err_reg_normal  = regionMeansLong(i1_reg_normal - i2_reg_normal, Constants.regions, 'Error', 'Registered');

i1_cata_normal_region_means = regionMeansLong(i1_cata_normal, Constants.regions, 'I1', 'Cata');
i1_1mid_normal_region_means = regionMeansLong(i1_1mid_normal, Constants.regions, 'I1', '1 Midlines');
i1_2mid_normal_region_means = regionMeansLong(i1_2mid_normal, Constants.regions, 'I1', '2 Midlines');
i1_reg_normal_region_means = regionMeansLong(i1_reg_normal, Constants.regions, 'I1', 'Registered');

i2_cata_normal_region_means = regionMeansLong(i2_cata_normal, Constants.regions, 'I2', 'Cata');
i2_1mid_normal_region_means = regionMeansLong(i2_1mid_normal, Constants.regions, 'I2', '1 Midlines');
i2_2mid_normal_region_means = regionMeansLong(i2_2mid_normal, Constants.regions, 'I2', '2 Midlines');
i2_reg_normal_region_means = regionMeansLong(i2_reg_normal, Constants.regions, 'I2', 'Registered');

% dist_means_dupes = regionMeansLong(d_normal, Constants.regions, 'Distance', 'NA');

errors_table = vertcat(err_cata_normal, err_1mid_normal, err_2mid_normal, err_reg_normal);
I1_table = vertcat(i1_cata_normal_region_means, i1_1mid_normal_region_means, i1_2mid_normal_region_means, i1_reg_normal_region_means);
I2_table = vertcat(i2_cata_normal_region_means, i2_1mid_normal_region_means, i2_2mid_normal_region_means, i2_reg_normal_region_means);

big_table = join(errors_table, join(I1_table, I2_table));

% big_table.Distance = repmat(dist_means_dupes.Distance, [4 1]);
% % big_table = join(dist_means_normal, join(err_1mid_normal, join(err_2mid_normal, join(err_reg_normal, err_cata_normal))));
writetable(big_table, '~/Desktop/normal_errors_table.csv');
%% Plot Normal
% data_ = cat(3, ...
%     i1_cata_normal ./ i2_cata_normal,...
%     i1_1mid_normal ./ i2_1mid_normal,...
%     i1_2mid_normal ./ i2_2mid_normal,...
%     i1_reg_normal ./ i2_reg_normal);


data_ = cat(3,...
    abs(i1_cata_normal - i2_cata_normal),...
    abs(i1_1mid_normal - i2_1mid_normal),...
    abs(i1_2mid_normal - i2_2mid_normal),...
    abs(i1_reg_normal - i2_reg_normal));

labels_ = {'Cata', '1 Midline, Two Masks (Unregistered)', '2 Midlines, Two Masks (Unregistered)', '2 Midlines, Two Masks (Registered)'};
ylim_ = [0 1000];
plotMultiplePharynxData(data_, labels_, ylim_);
addRegionBoundsToPlot(gca, Constants.regions);

%%
data_ = cat(3,...
    abs(i1_cata_normal - i2_cata_normal),...
    abs(i1_1mid_normal - i2_1mid_normal),...
    abs(i1_2mid_normal - i2_2mid_normal),...
    abs(i1_reg_normal - i2_reg_normal));

labels_ = {'Cata', '1 Midline, Two Masks (Unregistered)', '2 Midlines, Two Masks (Unregistered)', '2 Midlines, Two Masks (Registered)'};
ylim_ = [0 1000];
plotMultiplePharynxData(data_, labels_, ylim_);
addRegionBoundsToPlot(gca, Constants.regions);
%% PLOTS
data_ = cat(3, ...
    i1_cata_normal,...
    i1_1mid_normal,...
    i1_2mid_normal,...
    i1_reg_normal);
labels_ = {'Cata', '1 Midline, Two Masks (Unregistered)', '2 Midlines, Two Masks (Unregistered)', '2 Midlines, Two Masks (Registered)'};
ylim_ = [0 13000];
plotMultiplePharynxData(data_, labels_, ylim_);
addRegionBoundsToPlot(gca, Constants.regions);

%% Helpers
function pairs = unorderedPairs(v)
    [A,B] = meshgrid(v,v);
    c=cat(2,A',B');
    pairs=reshape(c,[],2);
%     pairs(pairs(:,1)==pairs(:,2),:) = [];
end