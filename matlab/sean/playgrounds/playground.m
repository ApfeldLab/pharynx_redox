%%
% Idea: generate fake bad data by taking the ratio of images at all other
% time points

pairs = combnk(1:nAnimals, 2);
fake_ratios = zeros(size(im410_animal1_frame1, 1), size(im410_animal1_frame1, 2), length(pairs));
fake_ratio_colored = zeros(size(im410_animal1_frame1, 1), size(im410_animal1_frame1, 2), 3, length(pairs)); % RGB
cmap = cbrewer('div', 'RdBu', 256, 'PCHIP');

min_ = 0.8681;
max_ = 1.1319;

for i=1:length(pairs)
    im1_idx = pairs(i,1);
    im2_idx = pairs(i,2);
    fake_ratios(:,:,i) = im410_animal1_frame1_raw(:,:,pairs(i,1)) ./ im410_animal1_frame1_raw(:,:,pairs(i,2));
    

    fake_ratio_colored(:,:,:,i) = ja_adjust_brightness(fake_ratios(:,:,i), im410_animal1_frame1(:,:,pairs(i,1)), 2500, cmap, min_, max_);
end
imshow(fake_ratio_colored(30:80,40:110,:,1))

%%

%% Write fake  ratio colored to GIF
S = fake_ratio_colored;
filename = '~/Desktop/fake_colored_worms.gif';
for f = 1:size(S,4)
  [SIf,cm] = rgb2ind(S(40:80,40:110,:,f),256);
  if f == 1
    imwrite(SIf,cm,filename,'Loop',Inf,'Delay',0);
  else
    imwrite(SIf,cm, filename,'WriteMode','append','Delay',0);
  end
end

%%
i = 2;
idata = i410_animal_1_frame_1(:,i);
fd_tmp = makeWormFd_SJ(idata, 'lambda', 10^-3);
fd_smooth = makeWormFd_SJ(idata);
fd_avg = makeWormFd_SJ(i410_animal_1_frame_1);
plot(idata, 'o'); hold on;
plot(fd_tmp);
plot(fd_smooth);
plot(fd_avg);

%%
sean_r = i410_animal_1_frame_1 ./ i410_animal_1_frame_2;
cata_r = i_animal_1_frame_1_cata ./ i_animal_1_frame_2_cata;
reg_r = eval_fd(1:100, reg_i410_animal1_frame1 ./ reg_i410_animal1_frame2);

hm = m_sep.AnteriorBulb == 3;

plotMultiplePharynxData(cat(3, cata_r(:,hm), sean_r(:,hm), reg_r(:,hm)), {'Sean', 'Cata', 'Reg'}, [.95 1.05]); 
hold on; hline(1);
addRegionBoundsToPlot(gca, Constants.regions);

%%
reg_i410_1 = eval_fd(1:100, reg_i410_animal1_frame1);
plotMultiplePharynxData(cat(3, i_animal_1_frame_1_cata(:,hm), i410_animal_1_frame_1(:,hm), reg_i410_1(:,hm)), {'Sean', 'Cata', 'Reg'}); 
addRegionBoundsToPlot(gca, Constants.regions);

%%
i_animal_1_frame_1_cata = measureIntensityCata(im410_animal1_frame1, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
i_animal_1_frame_2_cata = measureIntensityCata(im410_animal1_frame2, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
figure;
plot(i_animal_1_frame_1_cata, '-r'); hold on;
plot(i_animal_1_frame_2_cata, '-b');
%%
% i = 27;

figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:57
clf;
HIDPlot(...
    imR(:,:,i), seg410_animal1_frame1(:,:,i), ...
    im410_animal1_frame1(:,:,i), im410_animal1_frame2(:,:,i),...
    i410_animal_1_frame_1(:,i), i410_animal_1_frame_2_midline_1(:,i), i410_animal_1_frame_2(:,i),...
    i_animal_1_frame_1_cata(:,i), i_animal_1_frame_2_cata(:,i),...
    reg_i410_animal1_frame1(i), reg_i410_animal1_frame2(i), ...
    midlines410_animal1_frame1{i}, midlines410_animal1_frame2{i}, ...
    scaledLrBounds_animal_1_frame_1(i,:), scaledLrBounds_animal_1_frame_2(i,:),...
    fdObjs_animal1(i).warpFD, m_sep(i, :), i);
export_fig(sprintf("/Users/sean/Desktop/HID Plots/%d.pdf", i));
end
%%
single_mask_i410_1 = measureIntensityCata(im410_animal1_frame1, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
single_mask_i410_2 = measureIntensityCata(im410_animal1_frame2, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);

double_mask_i410_1 = measureIntensityCata(im410_animal1_frame1, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
double_mask_i410_2 = measureIntensityCata(im410_animal1_frame1, seg410_animal1_frame2, midlines410_animal1_frame1, 1000);

cata_iR = measureIntensityCata(imR, seg410_animal1_frame1, midlines410_animal1_frame1, 1000);
%%
i=52;
figure;
hold on;
plot(single_mask_i410_1(:,i));
plot(single_mask_i410_2(:,i));
yyaxis right;
plot(cata_iR(:,i), 'r-');
plot(single_mask_i410_1(:,i) ./ single_mask_i410_2(:,i), 'b-');
hline(1);
ylim([.8 1.2]);
legend('I_{410_1}', 'I_{410_2}', 'I_R', 'I_{410_1} / I_{410_2}');
title('Single Mask');
hold off;

figure;
hold on;
plot(double_mask_i410_1(:,i));
plot(double_mask_i410_2(:,i));
yyaxis right;
plot(cata_iR(:,i), 'r-');
plot(double_mask_i410_1(:,i) ./ double_mask_i410_2(:,i), 'b-');
hline(1);
ylim([.8 1.2]);
legend('I_{410_1}', 'I_{410_2}', 'I_R', 'I_{410_1} / I_{410_2}');
title('Double Mask');
hold off;

%%
figure;
ax = gca;
% for a=1:size(im410_animal1_frame1, 3)
    a=52;
    plotMultipleMidlines(imR(:,:,a), seg410_animal1_frame1(:,:,a), midlines410_animal1_frame1{a}, midlines410_animal1_frame2{a}, scaledLrBounds_animal_1_frame_1(a,:),scaledLrBounds_animal_1_frame_2(a,:), fdObjs_animal1(a).warpFD, ax);
    title(ax, num2str(a));
%     waitforbuttonpress;
% end

%%
imshow(cat(3,im410_animal1_frame1(:,:,52),im410_animal1_frame2(:,:,52),zeros(size(im410_animal1_frame1(:,:,52)))));
%%
figure;
subplot(3,1,1);
imshow(ja_adjust_brightness(imR(:,:,a), max(im410_animal1_frame2(:,:,a), im410_animal1_frame1(:,:,a)), 1000, cmap, .85, 1.15));  hold on;plot(midlines410_animal1_frame1{a}); plot(midlines410_animal1_frame2{a});

subplot(3,1,2);
imshow(im410_animal1_frame1(:,:,a)); hold on; plot(midlines410_animal1_frame1{a});

subplot(3,1,3);
imshow(ja_adjust_brightness(imR(:,:,a), max(im410_animal1_frame2(:,:,a), im410_animal1_frame1(:,:,a)), 1000, cmap, .85, 1.15));  hold on;plot(midlines410_animal1_frame1{a}); plot(midlines410_animal1_frame2{a});


%%

figure;
plot(i410_animal_1_frame_1(:,27)); 
hold on;
plot(i410_animal_1_frame_2(:,27));

yyaxis right;
plot(linspace(1,100,1000), cata_iR(:,i), 'r-');
plot(i410_animal_1_frame_1(:,27) ./ (i410_animal_1_frame_2(:,27)), 'b-');

plot(reg_i410_animal1_frame1(27) ./ reg_i410_animal1_frame2(27));


legend('I_{410_1}', 'I_{410_2}', 'I_R', 'I_{410_1} / I_{410_2}', 'R_1 / R_2');
hline(1);
ylim([.8 1.2]);



hold off;
%% Load Data
rawImgFilePath = "/Users/sean/Desktop/meeting_12_6_18/data/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
allImgs = loadTiffImageStack(rawImgFilePath);

I = loadIndexer('~/Desktop/indexer_example.csv');
m = csvread("~/Desktop/meeting_12_6_18/data/movement.csv");
m = m(1:57);

m_sep = readtable("~/Desktop/meeting_12_6_18/data/movement_separated.csv");

%
im410_1_raw = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 1 & I.Strain == "Animal 1"));
im410_2_raw = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "410/410" & I.SetFrame == 2 & I.Strain == "Animal 1"));

nAnimals = size(im410_1_raw, 3);
imTL = allImgs(:,:,I.ImgFrame(I.LambdaGroup == "TL" & I.Strain == "Animal 1"));
imTL = repmat(imTL, [1 1 size(im410_1_raw, 3)]);

im410_1 = subtractMedians(im410_1_raw);
im410_2 = subtractMedians(im410_2_raw);

%
imR_animal_1 = im410_1_raw ./ im410_2_raw;

outputFileName = '/Users/sean/Desktop/meeting_12_6_18/data/R_animal_1.tif';
for K=1:length(imR_animal_1(1, 1, :))
   I = imR_animal_1(:,:,K);
   scaledI = (I-min(I(:))) ./ (max(I(:)-min(I(:))));
   imwrite(scaledI, outputFileName, 'WriteMode', 'append');
end

%
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

% Register
[ri410_1_both, ri410_2_both, ~, ~, ~] = ...
    ChannelRegister(i410_1_mid_410_1, i410_2_mid_410_2, 100);
[ri410_1_1st, ri410_2_1st, ~, ~, ~] = ...
    ChannelRegister(i410_1_mid_410_1, i410_2_mid_410_1, 100);
[ri410_1_2nd, ri410_2_2nd, ~, ~, ~] = ...
    ChannelRegister(i410_1_mid_410_2, i410_2_mid_410_2, 100);

%%
r_1st = i410_1_mid_410_1 ./ i410_2_mid_410_1;
r_2nd = i410_1_mid_410_2 ./ i410_2_mid_410_2;
r_both = eval_fd(1:100, ri410_1_both ./ ri410_2_both);

%% R midline comparison
x=1:100;
y_all = mean(r_both, 2);
y_1st = mean(r_1st, 2);
y_2nd = mean(r_2nd, 2);

SEM_all = std(r_both, 0, 2)/sqrt(length(y_all));
SEM_1st = std(r_1st, 0, 2)/sqrt(length(r_1st));
SEM_2nd = std(r_2nd, 0, 2)/sqrt(length(r_2nd));

ts = tinv([0.025  0.975],length(y_all)-1);
CI_both = y_all + ts.*SEM_all;
CI_1st = y_1st + ts.*SEM_1st;
CI_2nd = y_2nd + ts.*SEM_2nd;

e_both = abs(y_all-CI_both);
e_1st = abs(y_1st-CI_1st);
e_2nd= abs(y_2nd-CI_2nd);
figure;
set(gcf, 'Position', get(0, 'Screensize'));
cmap_scaled = flipud(cbrewer('qual', 'Dark2', 3));
[l, p] = boundedline(x, y_1st, e_1st, ' ', ...
    x, y_2nd, e_2nd, ...,
    x, y_all, e_both, 'cmap', cmap_scaled);
outlinebounds(l, p);
xlim([1 100]);
ylim([.98 1.02]);
legend('1st/1st (Trimmed)', '2nd/2nd (Trimmed)', '1st/2nd (Registered)', 'Location', 'northeast');
title('Effect of Midline Strategy on R_{410/410}');
xlabel('Pharynx Space');
ylabel('R_{410/410}');

export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_strategies.pdf
close all;

%% R Movement Comparison, Both Midlines
y_all = mean(r_both, 2);
SEM_all = std(r_both, 0, 2)/sqrt(length(y_all));
ts = tinv([0.025  0.975],length(y_all)-1);
CI_both = y_all + ts.*SEM_all;

y_others = zeros(5, length(y_all));

for mvmt_code=0:4
    
end

%% Register
[reg410_FD, reg470_FD, warp470, resampled_intensity, fdObjs] = ...
    ChannelRegister(i410_1_mid_410_1, i410_2_mid_410_2, 100);

%% Show image with colormap
cmap = cbrewer('div', 'RdBu', 1000);
%%
cmap_scaled = MC_nonlinearCmap(cmap, 1, [.5 1.5], 1, 1);
%% `Info` Plot
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
    plot(r_1st_midline_only(:,i));
    plot(r_2nd_midline_only(:,i));
    plot(r_both(:,i), 'Color', 'black');
    ylim([low hi]);
    legend('1st/1st', '2nd/2nd', '1st/2nd');
    title('R_{410/410}(s) for Midline Strategies');
    
    % I(s)
    subplot(2,3,4:5);
    hold on;
    plot(i410_1_mid_410_1(:,i));
    plot(i410_2_mid_410_2(:,i));
    plot(1:100, eval_fd(1:100, ri410_2_both(i)), 'Color', 'black');
    ylim([0 15000]);
    legend('1st', '2nd', '2nd (reg)');
    title('I_{410}(s)');
    
    % R(x,y)
    subplot(2,3,3);
    I = imR_animal_1(:,:,i);
    I_masked = I.*seg410_1(:,:,i);
    I_masked(I_masked == 0) = NaN;
    imAlpha=ones(size(I_masked));
    imAlpha(isnan(I_masked))=0;

    bbox = regionprops(seg410_1(:,:,1), 'BoundingBox');
    bbox = floor(bbox.BoundingBox);

    buffer = 15;

    left = bbox(1) - buffer;
    top = bbox(2) - buffer;
    right = bbox(1) + bbox(3) + buffer;
    bottom = bbox(2) + bbox(4) + buffer;
    imagesc(I_masked(top:bottom,left:right), 'AlphaData', imAlpha(top:bottom,left:right)); 
    colormap(cmap); colorbar;
    caxis([low hi]);
    set(gca,'color',0*[1 1 1]);
    axis square;
    title('R_{410/410}(x,y)');
    
    % Overall Figure
    sgtitle(sprintf('Animal %d (%d)', i, m(i)));
    
    % Save
%     export_fig('/Users/sean/Desktop/meeting_12_6_18/figures/ratio_images/info.pdf', '-opengl', '-append');
    waitforbuttonpress;
end

%%
figure;
ax=gca;
set(ax,'Color','k');
% for i=1:size(seg_410, 3)
for i=[4, 5, 27, 33]
    plotMultipleMidlines(im_R(:,:,i), seg_410(:,:,i), midlines_410{i}, midlines_470{i}, scaledLRBounds_410(i,:), scaledLRBounds_470(i,:), fdObjs(i).warpFD, ax);
    export_fig(sprintf("/Users/sean/Desktop/midline_plots_painters/%d.pdf", i), "-painters");
    cla(ax);
end

%% Movement Scatter Plot
cmap = cbrewer('qual', 'Set1', 3);
Y = table2array(m_sep);
jitterAmount = 0.07;
jitterValuesX = 2*(rand(size(Y(:,1)))-0.5)*jitterAmount;   % +/-jitterAmount max
jitterValuesY = 2*(rand(size(Y(:,2)))-0.5)*jitterAmount;   % +/-jitterAmount max
jitterValuesZ = 2*(rand(size(Y(:,3)))-0.5)*jitterAmount;   % +/-jitterAmount max

figure;
scatter3(Y(:,1) + jitterValuesX, Y(:,2) + jitterValuesY, Y(:,4) + jitterValuesZ, 50, m>0);
colormap(cmap);
colorbar;
xlabel('Posterior');
ylabel('Anterior');
zlabel('Tip');
title('New Movment Classification (Colored by Old Classification)');

%%

% Both
R = r_both;
plotRForMovement(R, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (Both Midlines)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/both/R_posterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (Both Midlines)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/both/R_anterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.SidesOfTip, 'R for Side Movement (Both Midlines)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/both/R_side_movement.pdf;
close all;
plotRForMovement(R, m_sep.Tip, 'R for Tip Movement (Both Midlines)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/both/R_tip_movement.pdf;
close all;

R = r_1st_midline_only;
plotRForMovement(R, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (1st Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/first/R_posterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (1st Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/first/R_anterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.SidesOfTip, 'R for Side Movement (1st Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/first/R_side_movement.pdf;
close all;
plotRForMovement(R, m_sep.Tip, 'R for Tip Movement (1st Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/first/R_tip_movement.pdf;
close all;

R = r_2nd_midline_only;
plotRForMovement(R, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (2nd Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/second/R_posterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (2nd Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/second/R_anterior_movement.pdf;
close all;
plotRForMovement(R, m_sep.SidesOfTip, 'R for Side Movement (2nd Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/second/R_side_movement.pdf;
close all;
plotRForMovement(R, m_sep.Tip, 'R for Tip Movement (2nd Midline Only)');
export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_per_region_movement/second/R_tip_movement.pdf;
close all;

%% 
ax = subplot(4, 3, 1);
plotRForMovement(r_both, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (1st/2nd)', ax);
ax = subplot(4, 3, 2);
plotRForMovement(r_1st, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (1st/1st)', ax);
ax = subplot(4, 3, 3);
plotRForMovement(r_2nd, m_sep.PosteriorBulb, 'R for Posterior Bulb Movement (2nd/2nd)', ax);
legend('0', '1', '2', '3', 'Location', 'southeast', 'Orientation', 'horizontal');

ax = subplot(4, 3, 4);
plotRForMovement(r_both, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (1st/2nd)', ax);
ax = subplot(4, 3, 5);
plotRForMovement(r_1st, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (1st/1st)', ax);
ax = subplot(4, 3, 6);
plotRForMovement(r_2nd, m_sep.AnteriorBulb, 'R for Anterior Bulb Movement (2nd/2nd)', ax);

ax = subplot(4, 3, 7);
plotRForMovement(r_both, m_sep.AnteriorBulb, 'R for Side Movement (1st/2nd)', ax);
ax = subplot(4, 3, 8);
plotRForMovement(r_1st, m_sep.AnteriorBulb, 'R for Side Movement (1st/1st)', ax);
ax = subplot(4, 3, 9);
plotRForMovement(r_2nd, m_sep.AnteriorBulb, 'R for Side . Movement (2nd/2nd)', ax);

ax = subplot(4, 3, 10);
plotRForMovement(r_both, m_sep.AnteriorBulb, 'R for Tip Movement (1st/2nd)', ax);
ax = subplot(4, 3, 11);
plotRForMovement(r_1st, m_sep.AnteriorBulb, 'R for Tip Movement (1st/1st)', ax);
ax = subplot(4, 3, 12);
plotRForMovement(r_2nd, m_sep.AnteriorBulb, 'R for Tip Movement (2nd/2nd)', ax);

%%
%%
jitterAmount = 1;
jitterValuesX = 2*(rand(size(Y(:,1)))-0.5)*jitterAmount;   % +/-jitterAmount max
jitterValuesY = 2*(rand(size(Y(:,2)))-0.5)*jitterAmount;   % +/-jitterAmount max

c = getcoef(ri410_1_both).';
Y = tsne(c);
figure;gscatter(Y(:,1)+jitterValuesX, Y(:,2)+jitterValuesY,m);

%% 
function plotRForMovement(r, mvmt_mode, t, ax)
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

% figure;
% set(gcf, 'Position', get(0, 'Screensize'));
cmap_scaled = flipud(cbrewer('qual', 'Dark2', 4));
[l, p] = boundedline(x, y_data, e_data, 'cmap', cmap_scaled, ax);
outlinebounds(l, p);
xlim([1 100]);
ylim([.95 1.05]);
% legend('0', '1', '2', '3', 'Location', 'northeast');
title(t);
xlabel('Pharynx Space');
ylabel('R_{410/410}');

% export_fig /Users/sean/Desktop/meeting_12_6_18/figures/r_strategies.pdf
% close all;
end

%%



%% HELPER FUNCTIONS
function [i, unscaled_bounds, scaled_bounds] = measureAndTrim(imStack, midlines, lrBounds, profileLength)
    i_raw = measureIntensityAlongMidlines(imStack, midlines, lrBounds, profileLength, 'BILINEAR');
    [i_trimmed, unscaled_bounds, scaled_bounds] = trimProfile(i_raw, lrBounds);
    i = ssquare(i_trimmed);
end