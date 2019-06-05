imageStackFileName = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif";

experimentDirectory = fileparts(imageStackFileName);
I = loadStandardIndexer(fullfile(experimentDirectory, "indexer.csv"));
storageDirectory = fullfile(experimentDirectory, "analysis/");

%% Load Movement Data
movement_pair1 = readtable(fullfile(experimentDirectory, "aggregated_movement_annotations_pair1.csv"));
movement_pair2 = readtable(fullfile(experimentDirectory, "aggregated_movement_annotations_pair2.csv"));

movement_pair1.Properties.VariableNames = {'frame1', 'posterior1', 'anterior1', 'sides_of_tip1', 'tip1'};
movement_pair2.Properties.VariableNames = {'frame2', 'posterior2', 'anterior2', 'sides_of_tip2', 'tip2'};

I = join(I, movement_pair1, 'LeftKeys', {'Animal'}, 'RightKeys', {'frame1'});
I = join(I, movement_pair2, 'LeftKeys', {'Animal'}, 'RightKeys', {'frame2'});

%%
allImages = loadTiffImageStack(imageStackFileName);
imTL    = allImages(:,:,I.ImgFrame(I.Lambda == "TL"));
im410_1 = allImages(:,:,I.ImgFrame(I.Lambda == "410_1"));
im470_1 = allImages(:,:,I.ImgFrame(I.Lambda == "470_1"));
im410_2 = allImages(:,:,I.ImgFrame(I.Lambda == "410_2"));
im470_2 = allImages(:,:,I.ImgFrame(I.Lambda == "470_2"));

[i1_pair1, i2_pair1] = pipelineTwoMidlinesTwoMasks(imTL, im410_1, im470_1);
[i1_pair2, i2_pair2] = pipelineTwoMidlinesTwoMasks(imTL, im410_2, im470_2);

r_pair1 = i1_pair1 ./ i2_pair1;
r_pair2 = i1_pair2 ./ i2_pair2;

OxD = ja_oxd(r_pair1);
E = ja_E(OxD);

%% Plots
strains = unique(I.Strain);

% Ratio by Strain
figure;
% ax1 = subplot(2,1,1);
plotMultiplePharynxData(separateDataByStrain(E, I), strains, [-280 -240], gca);

% ax2 = subplot(2,1,2);
% cmap = cbrewer('qual', 'Paired', 2*length(strains));
% [d, l] = separateDataByStrainAndMovement(r_pair1, I, 'posterior1');
% plotMultiplePharynxData(d, l, [1 2], ax2, cmap);

%% Scratchpad

centerx = size(rot_seg410, 2) / 2;
centery = size(rot_seg410, 1) / 2;
x = centerx - 15;
y = centery - 15;
w = 25;
h = 30;

imshow(fliplr(rot_seg410(:,:,1)), []);
rectangle('Position',[centerx - 15, centery - 15, 25, 30],...
          'EdgeColor', 'red',...
          'LineWidth',2);