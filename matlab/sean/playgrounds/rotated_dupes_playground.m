rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
allImgs = loadTiffImageStack(rawImgFilePath);
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
im410_1 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL"));

nAnimal1 = size(I.ImgFrame(I.LambdaGroup == "410/410" & I.Strain == "Animal 1" & I.SetFrame == 1), 1);
imTL = repmat(imTL(:,:,1), [1 1 nAnimal1]);

%% Rotate Images according to first frame
seg410_1 = segmentPharynx(im410_1, 0, 2000);

rot_410_1 = rotatePharynx(im410_1, seg410_1);
rot_410_2 = rotatePharynx(im410_2, seg410_1);
rot_TL = rotatePharynx(imTL, seg410_1);

%% Crop
crop_410_1 = rot_410_1(52:79, 54:121, :);
crop_410_2 = rot_410_2(52:79, 54:121, :);
crop_TL = rot_TL(52:79, 54:121, :);

%% Make New Image Stacks from All Permutation Pairs
pairsAnimal1 = unorderedPairs(1:nAnimal1);

im410_1_dupes = crop_410_1(:,:,pairsAnimal1(:,1));
im410_2_dupes = crop_410_2(:,:,pairsAnimal1(:,2));

imR_dupes = im410_1_dupes ./ im410_2_dupes;
imR = (crop_410_1 ./ crop_410_2);

save('~/Desktop/imR.mat', 'imR');
save('~/Desktop/imR_dupes.mat', 'imR_dupes');   

%% Helpers
function pairs = unorderedPairs(v)
    [A,B] = meshgrid(v,v);
    c=cat(2,A',B');
    pairs=reshape(c,[],2);
%     pairs(pairs(:,1)==pairs(:,2),:) = [];
end