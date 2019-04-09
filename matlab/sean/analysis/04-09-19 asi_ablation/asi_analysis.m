rawImgFilePath = "/Users/sean/code/wormAnalysis/data/2019_04_03_ASI_ablation/2019_04_03_ASI_ablation.tif";
allImgs = loadTiffImageStack(rawImgFilePath);
I = loadIndexer("/Users/sean/code/wormAnalysis/data/2019_04_03_ASI_ablation/indexer.csv");

%%
imTL = allImgs(:,:,I.ImgFrame(I.Lambda == "TL"));
im410 = allImgs(:,:,I.ImgFrame(I.Lambda == "410_1"));
im470 = allImgs(:,:,I.ImgFrame(I.Lambda == "470_1"));


% Registered
% [i1, i2, matchingVecs, midlines1, midlines2, scaled_bounds1, scaled_bounds2, ftdObjs, dx, dy, unreg_i1, unreg_i2] = pipelineTwoMidlinesTwoMasksRegistration(imTL, im410, im470);

% Unregistered
[i1, i2, seg1, seg2, midlines1, midlines2, lrBounds1, lrBounds2, i1_raw, i2_raw] = pipelineTwoMidlinesTwoMasks(imTL, im410, im470);

%
i_ratio = i1 ./ i2;

wtIDX = [1:20, 45:84];
asiIDX = [21:44, 85:123];

figure;
ax = subplot(2,2,1);
plotMultiplePharynxData({i1(:,wtIDX), i1(:,asiIDX)}, {"WT", "ASI"}, [0, 12000], ax);
title('I_{410}');

ax = subplot(2,2,2);
plotMultiplePharynxData({i2(:,wtIDX), i2(:,asiIDX)}, {"WT", "ASI"}, [0, 12000], ax);
title('I_{470}');

ax = subplot(2,2,3);
iR = i1./i2;
plotMultiplePharynxData({iR(:,wtIDX), iR(:,asiIDX)}, {"WT", "ASI"}, [.9, 1.35], ax);
title('I_{410/470}');

ax = subplot(2,2,4);
E = ja_E(ja_oxd(iR));
plotMultiplePharynxData({E(:,wtIDX), E(:,asiIDX)}, {"WT", "ASI"}, [-292, -272], ax);
title('E');
%%
figure;
plotMultipleMidlines(imR(:,:,1), im410(:,:,1), seg1(:,:,1), midlines1{1}, midlines2{1}, ...
    scaled_bounds1, scaled_bounds2, fdObjs(1).warpFD);