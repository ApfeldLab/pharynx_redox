%% Load Data
rawImgFilePath = "/Users/sean/Desktop/meeting_12_6_18/data/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
allImgs = loadTiffImageStack(rawImgFilePath);

I = loadIndexer('~/Desktop/indexer_example.csv');
m = csvread("~/Desktop/meeting_12_6_18/data/movement.csv");
m = m(1:57);

m_sep = readtable("~/Desktop/meeting_12_6_18/data/movement_separated.csv");


im410_1_raw = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2_raw = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));

nAnimals = size(im410_1_raw, 3);
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL" & I.Strain == "Animal 1"));
imTL = repmat(imTL, [1 1 size(im410_1_raw, 3)]);

im410_1 = subtractMedians(im410_1_raw);
im410_2 = subtractMedians(im410_2_raw);

seg410_1 = segmentPharynx(im410_1, 0);
seg410_2 = segmentPharynx(im410_2, 0);

segLRBounds_1 = double(getLeftRightBounds(seg410_1));
segLRBounds_2 = double(getLeftRightBounds(seg410_2));


midlines410_1 = calculateMidlines(imTL, seg410_1);
midlines410_2 = calculateMidlines(imTL, seg410_2);

% Trim
INTERP_METHOD = 'bilinear';
PROFILE_LENGTH = 1000;

[i410_1_mid_410_1, ~, scaledLrBounds_11] = measureAndTrim(im410_1, midlines410_1, segLRBounds_1, PROFILE_LENGTH, INTERP_METHOD);
[i410_1_mid_410_2, ~, scaledLrBounds_12] = measureAndTrim(im410_1, midlines410_2, segLRBounds_2, PROFILE_LENGTH, INTERP_METHOD);
[i410_2_mid_410_1, ~, scaledLrBounds_21] = measureAndTrim(im410_2, midlines410_1, segLRBounds_1, PROFILE_LENGTH, INTERP_METHOD);
[i410_2_mid_410_2, ~, scaledLrBounds_22] = measureAndTrim(im410_2, midlines410_2, segLRBounds_2, PROFILE_LENGTH, INTERP_METHOD);

% % Register
% [ri410_1_both, ri410_2_both, ~, ~, ~] = ...
%     ChannelRegister(i410_1_mid_410_1, i410_2_mid_410_2, 100);
% % [ri410_1_1st, ri410_2_1st, ~, ~, ~] = ...
% %     ChannelRegister(i410_1_mid_410_1, i410_2_mid_410_1, 100);
% % [ri410_1_2nd, ri410_2_2nd, ~, ~, ~] = ...
% %     ChannelRegister(i410_1_mid_410_2, i410_2_mid_410_2, 100);

%% Plot Errors

figure;
e_I = i410_1_mid_410_1 - i410_2_mid_410_1;
sum_e_I = sum(e_I);
sum_I = sum(i410_1_mid_410_1);
scatter(sum_I, sum_e_I);
ylabel('$\sum{\epsilon_I}$', 'Interpreter', 'latex');
xlabel('$\sum{I}$', 'Interpreter', 'latex');
title('The Error in Intensity is Correlated with Intensity');

figure;
sum_I1 = sum(i410_1_mid_410_1);
sum_I2 = sum(i410_2_mid_410_1);
scatter(sum_I1, sum_I2);
xlabel('$\sum{I_1}$', 'Interpreter', 'latex');
ylabel('$\sum{I_2}$', 'Interpreter', 'latex');
title('Cumulative Intensity is Correlated Across Frames');

%%
figure;
residual_e = sum_e_I ./ sum_I;
scatter(sum_I, residual_e);

%% Let's look only at pm7

% First, let's draw a plot of the average of the 1st intensity profile also
% adding the region boundaries
avg_I = mean(i410_1_mid_410_1, 2);
upper_lim = max(avg_I) * 1.1;
figure;
plot(avg_I);
ylim([0 upper_lim]);
% region boundaries
% Add a patch
fn = fieldnames(Constants.regions);
for i=1:numel(fn)-2
    b = Constants.regions.(fn{i});
    patch([b(1) b(1) b(2) b(2)],[0 upper_lim upper_lim 0], [0 0 0], 'FaceAlpha', 0.1);
    text(b(1)+.5, upper_lim - 500, fn{i});
end