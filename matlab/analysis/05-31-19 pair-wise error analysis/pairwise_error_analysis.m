    clear;clc;
imageFilePaths = {
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_23-HD233_HD236/2017_02_23-HD233_HD236.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_27-HD233_HD236/2017_02_27-HD233_HD236.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_03_01_HD233_HD236/2017_03_01_HD233_HD236.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_08_15-HD233_4mM-lev/2017_08_15-HD233_4mM-lev.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_08_23-HD233_4mm_lev/2017_08_23-HD233_4mm_lev.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_08_24-HD233_SAY93/2017_08_24-HD233_SAY93.tif"
    "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_08_25-HD233_4mm_lev/2017_08_25-HD233_4mm_lev.tif"
%     "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_08_25-HD233_SAY93/2017_08_25-HD233_SAY93.tif"
%     "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_10_12-gld_1_RNAi_SAY88_HD233/2017_10_12-gld_1_RNAi_SAY88_HD233.tif"
};

allImages   = cell(length(imageFilePaths), 1);
allIndexers = cell(length(imageFilePaths), 1);
experimentNames = cell(length(imageFilePaths), 1);

for i=1:length(imageFilePaths)
    imageFilePath = imageFilePaths{i};
    disp(imageFilePath);
    allImages{i} = loadTiffImageStack(imageFilePath);
    allIndexers{i} = loadStandardIndexer(fullfile(fileparts(imageFilePath), "indexer.csv"));
    
    [~, ename, ~] = fileparts(imageFilePath);
    experimentNames{i} = repelem(ename, height(allIndexers{i})).';
end

allImages = cat(3, allImages{:});
Indexer = vertcat(allIndexers{:});
Indexer.experiment = vertcat(experimentNames{:});
Indexer.Animal = ceil((1:height(Indexer))/5).';
Indexer.strain_experiment = strcat(Indexer.Strain, " ", Indexer.experiment);
%%
imTL    = allImages(:,:,1:5:end);
im470_1 = allImages(:,:,2:5:end);
im410_1 = allImages(:,:,3:5:end);
im470_2 = allImages(:,:,4:5:end);
im410_2 = allImages(:,:,5:5:end);

[i1_pair1, i2_pair1] = pipelineTwoMidlinesTwoMasks(imTL, im410_1, im470_1);
[i1_pair2, i2_pair2] = pipelineTwoMidlinesTwoMasks(imTL, im410_2, im470_2);

r_pair1 = i1_pair1 ./ i2_pair1;
r_pair2 = i1_pair2 ./ i2_pair2;

OxD = ja_oxd(r_pair2);
E = ja_E(OxD);

%% Plots

strains = unique(Indexer.Strain);
% strains = unique(Indexer.strain_experiment);

% Ratio by Strain
figure;
% ax1 = subplot(2,1,1);
cmap = cbrewer('qual', 'Set1', length(strains));
% plotMultiplePharynxData(separateDataByStrain(E, I), strains, [-320 -240], gca, cmap);

plotMultiplePharynxData(separateDataByStrain(fliplr(E), Indexer), strains, [-275 -265], gca, cmap);