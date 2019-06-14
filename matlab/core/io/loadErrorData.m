function [imTL,im410_1, im410_2,movement] = loadErrorData()
    %LOADERRORDATA Summary of this function goes here
    %   Detailed explanation goes here
    rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
    allImgs = loadTiffImageStack(rawImgFilePath);
    I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
    imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL"));
    im410_1 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
    im410_2 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));

    nAnimal1 = size(I.ImgFrame(I.LambdaGroup == "410/410" & I.Strain == "Animal 1" & I.SetFrame == 1), 1);
    imTL = repmat(imTL(:,:,1), [1 1 nAnimal1]);

    movement = readtable("~/Desktop/meeting_12_6_18/data/movement_separated.csv");
end

