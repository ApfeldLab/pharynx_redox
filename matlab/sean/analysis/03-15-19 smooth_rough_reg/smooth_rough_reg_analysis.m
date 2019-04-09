rawImgFilePath = "/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/WT high mvmt 11-19-18.tif";
allImgs = loadTiffImageStack(rawImgFilePath);
I = loadIndexer('/Users/sean/code/wormAnalysis/matlab/sean/analysis/01-21-19 error analysis/data/indexer.csv');
im410_1 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2 = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL"));

nAnimal1 = size(I.ImgFrame(I.LambdaGroup == "410/410" & I.Strain == "Animal 1" & I.SetFrame == 1), 1);
imTL = repmat(imTL(:,:,1), [1 1 nAnimal1]);

% Rotate Images according to first frame
seg410_1 = segmentPharynx(im410_1, 0, 2000);

im410_1 = rotatePharynx(im410_1, seg410_1);
im410_2 = rotatePharynx(im410_2, seg410_1);
imTL = rotatePharynx(imTL, seg410_1);

% Crop
% im410_1 = im410_1(52:79, 54:121, :);
% im410_2 = im410_2(52:79, 54:121, :);
% imTL = imTL(52:79, 54:121, :);

%
[i1, i2, matchingVecs, midlines1, midlines2, scaled_bounds1, scaled_bounds2, fdObjs, dx, dy, unreg_i1, unreg_i2] = pipelineTwoMidlinesTwoMasksRegistration(imTL, im410_1, im410_2);

%% Hyperparameter search
% Hyperparameters:
% Smooth Lambda, Rough Lambda, Warp Lambda

smooth_lambda_range = linspace(1e-1, 1e1, 10);
warp_lambda_range = linspace(1e-1, 1e1, 10);

combs = allcomb(smooth_lambda_range, warp_lambda_range);

smooth_lambda_combs = combs(:,1);
warp_lambda_combs   = combs(:,2);

mmaes = zeros(size(combs,1), 1);
xs = linspace(1,100,1000);
textprogressbar('Parameter Sweep: ');
for j=1:size(combs,1),
    textprogressbar(j);
    [~, ~, roughFD1_, ~, regRoughFD2_, ~, ~] = smoothRoughRegister(unreg_i1, unreg_i2, smooth_lambda_combs(j), 1e-3, warp_lambda_combs(j));
    ys_1 =eval_fd(xs, roughFD1_);
    ys_2 = eval_fd(xs, regRoughFD2_);
    mmaes(j) = mean(sum(abs(ys_1 - ys_2),2)/size(ys_1,1));
end
textprogressbar('Done');

%%
% figure;
% hold on;
% cdfplot(sum(abs(eval_fd(1:100, roughFD1) - eval_fd(1:100, roughFD2)),2)/size(ys_1,1)); hold on;
% cdfplot(sum(abs(eval_fd(1:100, roughFD1) - eval_fd(1:100, roughRegF2)),2)/size(ys_1,1)); hold on;
% hold off;
% legend('unregistered', 'registered');

