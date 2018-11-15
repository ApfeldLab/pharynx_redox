%% Load Images
rawImgFilePath = "/Users/sean/Documents/worm stuff/vab1_2do_05_25_mvmt_new/vab1_2do_05_25_mvmt.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
imTable = struct2table(tiffread2(rawImgFilePath));
allImgs = cat(3, imTable.data{:}); % dimensions: WxHxi; e.g. first image is ims(:, :, 1)
    
% TODO - Un-hardcode three channels
imTL = allImgs(:, :, 1:3:end);
im410_raw = double(allImgs(:, :, 2:3:end));
im470_raw = double(allImgs(:, :, 3:3:end));
nAnimals = size(im410_raw, 3);

im410_raw(im410_raw == 0) = 1;
im470_raw(im470_raw == 0) = 1;

im_r_410_470 = im410_raw ./ im470_raw;

% Load Movement
mvmnt = csvread(fullfile(rawImgDir, 'movement.csv'));

% Subtract medians
im410 = subtractMedians(im410_raw);
im470 = subtractMedians(im470_raw);

%% Strain Names
strainTable = readtable(fullfile(rawImgDir, 'strains.csv'));
strains = cell(1,nAnimals);
for i=1:size(strainTable, 1)
    row_ = strainTable(i,:);
    s = str2double(cell2mat(strainTable(i,:).Start));
    e = str2double(cell2mat(strainTable(i,:).End));
    strains(s:e) = row_.Strain;
end
strains = categorical(cellstr(strains.'));
[GN, ~, G] = unique(strains);

%% Process Images

% This segmentation usually works, if the data looks a little wonky,
% take a look at the segmentation. That's usually the problem. If the
% segmentation fails. Go to ImageJ and try segmenting there, then load the
% masks here instead of calling segmentPharynx.
seg410 = segmentPharynx(im410, 0);
seg470 = segmentPharynx(im470, 0);

[centroids, orientations] = calcCentroidsAndOrientations(seg410);

rotatedTL = rotatePharynx(imTL, centroids, orientations);
rotated410 = rotatePharynx(im410, centroids, orientations);
rotated470 = rotatePharynx(im470, centroids, orientations);

rotatedSeg410 = logical(rotatePharynx(seg410, centroids, orientations));
rotatedSeg470 = logical(rotatePharynx(seg470, centroids, orientations));

%% Bounding Boxes
lrBounds = getLeftRightBounds(seg410);

%% Generate Midlines
midlines = calculateMidlines(imTL, seg410);
midlines_rot = calculateMidlines(rotatedTL, rotatedSeg410);

%% get profile data
INTERP_METHOD = 'bilinear';
PROFILE_LENGTH = 1000;
i410 = measureIntensityAlongMidlines(im410, midlines, PROFILE_LENGTH, INTERP_METHOD);
i470 = measureIntensityAlongMidlines(im470, midlines, PROFILE_LENGTH, INTERP_METHOD);
i_r_410_470 = measureIntensityAlongMidlines(im_r_410_470, midlines, PROFILE_LENGTH, INTERP_METHOD);

i410_trimmed = ssquare(trimProfile(i410));
i470_trimmed = ssquare(trimProfile(i470));

%% Plot Strain Means

% Make a figures directory
mkdir(rawImgDir, 'figures');

plotStrainMeans(strains, i410_trimmed, "Mean 410nm \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_410nm_raw.png'));
plotStrainMeans(strains, i470_trimmed, "Mean 470nm \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_470nm_raw.png'));

% Calc R, OxD, E
R = i410_trimmed ./ i470_trimmed;
OxD = ja_oxd(R);
E = ja_E(OxD);

plotStrainMeans(strains, E, "Mean E \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_E_raw.png'));
plotStrainMeans(strains, OxD, "Mean OxD \pm std per Strain");
export_fig(fullfile(rawImgDir, 'figures', 'mean_OxD_raw.png'));

%% Save Profile Data
csvwrite(fullfile(rawImgDir, 'i410.csv'),i410_trimmed.');
csvwrite(fullfile(rawImgDir, 'i470.csv'),i470_trimmed.');

%% Plot according to Movment

% For example, filter out all >0s, then plot according to strain
E_no_movement = E(mvmnt==0);
strains_no_movement = strains(mvmnt==0);
plotStrainMeans(strains_no_movement, E, "");

%% Plot Image with Midline

i = 90;
I = rotated410(:,:,i);
midline = midlines_rot{i};
imshow(I, []); hold on;
plot(midline);
% rectangle('Position', bboxes(i, :),'EdgeColor','green');
hold off;

%% Plot normal line of midline at fixed x
i = 90;
I = masked_im410(:,:,i);
midline = midlines{i};

nLines = 50;
xs = linspace(double(lrBounds(i,1)), double(lrBounds(i,2)), nLines);
d_midline = differentiate(midline, xs);
ys_midline = feval(midline, xs);

imshow(I, []); hold on;
plot(midline);
for j=1:size(xs,2)
    x_j = xs(j);
    n_xs = linspace(x_j - 10, x_j + 10, 3); % These are the x-values of the normal line at x_j, (+/- 10)
    
    n_ys = (-1/d_midline(j)) .* (n_xs - x_j) + ys_midline(j);
    plot(n_xs, n_ys);
end

hold off;

%% measure means of normal lines

masked_ratio = seg410 .* im_r_410_470;
masked_ratio(masked_ratio == 0) = NaN;

masked_im410 = seg410 .* im410;
masked_im410(uint32(masked_ratio) == 0) = NaN;

nPoints = 1000;
xs = linspace(double(lrBounds(i,1)), double(lrBounds(i,2)), nPoints);

normal_means = zeros(nAnimals, nPoints);
normal_stds = zeros(nAnimals, nPoints);

for i=1:nAnimals
    I = masked_im410(:,:,i);
    midline = midlines{i};
    xs = linspace(double(lrBounds(i,1)), double(lrBounds(i,2)), nPoints);
    d_midline = differentiate(midline, xs);
    ys_midline = feval(midline, xs);
    for j=1:nPoints
        x_j = xs(j);
        n_xs = linspace(x_j - 10, x_j + 10, 20); % These are the x-values of the normal line at x_j, (+/- 10)
        n_ys = (-1/d_midline(j)) .* (n_xs - x_j) + ys_midline(j);

        prof_ = improfile(I, n_xs, n_ys, 10, 'bilinear');
        normal_means(i, j) = nanmean(prof_(:));
        normal_stds(i, j) = nanstd(prof_(:));
    end
end
normal_means = ssquare(trimProfile(normal_means.'));
%%
plot(normal_means(:,1)); hold on; plot(ssquare(i410_trimmed(:,1)));
legend('normal mean', 'midline');

%% measure variances of normal lines in Ratio Image for Movement Detection

masked_im410 = seg410 .* im410;
masked_im410(uint32(masked_im410) == 0) = NaN;

masked_ratio = seg410 .* im_r_410_470;
masked_ratio(uint32(masked_im410) == 0) = NaN;

nPoints = 20;
xs = linspace(double(lrBounds(i,1)), double(lrBounds(i,2)), nPoints);

normal_means = zeros(nAnimals, nPoints);
normal_stds = zeros(nAnimals, nPoints);
normal_vars = zeros(nAnimals, nPoints);

for i=1:nAnimals
    I = masked_ratio(:,:,i);
    midline = midlines{i};
    xs = linspace(double(lrBounds(i,1)), double(lrBounds(i,2)), nPoints);
    d_midline = differentiate(midline, xs);
    ys_midline = feval(midline, xs);
    for j=1:nPoints
        x_j = xs(j);
        n_xs = linspace(x_j - 10, x_j + 10, 20); % These are the x-values of the normal line at x_j, (+/- 10)
        n_ys = (-1/d_midline(j)) .* (n_xs - x_j) + ys_midline(j);

        prof_ = improfile(I, n_xs, n_ys, 10, 'bilinear');
        normal_means(i, j) = nanmean(prof_(:));
        normal_stds(i, j) = nanstd(prof_(:));
        normal_vars(i, j) = nanvar(prof_(:));
    end
end
normal_means_labels = horzcat([normal_means mvmnt]);
normal_vars_labels = horzcat([normal_stds mvmnt]);
normal_stds_labels = horzcat([normal_stds mvmnt]);