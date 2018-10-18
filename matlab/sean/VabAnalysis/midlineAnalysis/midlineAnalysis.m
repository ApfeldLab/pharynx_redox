%% Load data
DIRECTORY = '/Users/sean/Desktop/vab1_2do_05_25_mvmt';
COORDS_410_FNAME = 'PA-vab1_2do_05_25_mvmt_410_coords.dat';
COORDS_470_FNAME = 'PA-vab1_2do_05_25_mvmt_470_coords.dat';

INT_410_FNAME = 'PA-vab1_2do_05_25_mvmt_410_intensities.txt';
INT_470_FNAME = 'PA-vab1_2do_05_25_mvmt_470_intensities.txt';

coords410 = loadCoordinates(fullfile(DIRECTORY, COORDS_410_FNAME));
coords470 = loadCoordinates(fullfile(DIRECTORY, COORDS_470_FNAME));

mvmnt_rep = readtable(fullfile(DIRECTORY, 'movement.csv'), 'Format', '%d%d%*s');
mvmnt_rep = table2array(mvmnt_rep(:,2));

% Load images
r470410_metadata = tiffread2('/Users/sean/Desktop/vab1_2do_05_25_mvmt/vab1_2do_05_25_mvmt_R_470_410.tif');
r470410 = zeros(40, 90, size(r470410_metadata, 2));
for i = 1:size(r470410_metadata, 2)
    r470410(:, :, i) = r470410_metadata(i).data;
end

bin410_metadata = tiffread2('/Users/sean/Desktop/vab1_2do_05_25_mvmt/Mask of PA-vab1_2do_05_25_mvmt_410-1.tif');
bin410 = zeros(40, 90, size(bin410_metadata, 2));
for i = 1:size(bin410_metadata, 2)
    bin410(:, :, i) = bin410_metadata(i).data;
end
%%
nAnimals = size(coords410.x, 2);

%%
plotMidlines = 1;
f = figure;
set(f, 'Position', [90 301 1120 420], 'visible', 'off')
ax = axes('NextPlot', 'add');

xlim(ax, [5 90]);
ylim(ax, [1 40]);
textprogressbar('generating plots: ');
for animal = 1:nAnimals
    textprogressbar(i);
% animal = 1;
    cla(ax);
    imagesc(ax, r470410(:, :, animal));
    ax.Title.String = sprintf('Concomitant Points in the Channel-Specific Midlines (Animal %d)', animal);
    
    text(ax, 10, 14, sprintf('Movment: %d', mvmnt_rep(animal)), 'FontSize', 25);

    % Plot midlines
    if plotMidlines
        plot(ax, ...
            coords410.x(:,animal), coords410.y(:,animal), 'b', ...
            coords470.x(:,animal), coords470.y(:,animal), 'r');
    end
    % Plot arrows
    componentsU = coords470.x(:,animal) - coords410.x(:,animal);
    componentsV = coords470.y(:,animal) - coords410.y(:,animal);
    quiver(coords410.x(:,animal), coords410.y(:, animal), componentsU, componentsV, 'MaxHeadSize', 0.01, 'AutoScale', 0, 'Color', 'black');
    
    % Set legend
    if plotMidlines
        legend(ax, '410nm \rightarrow 470nm', '410nm', '470nm');
    else
        legend(ax, '410nm \rightarrow 470nm');
    end
    export_fig(sprintf('/Users/sean/Desktop/wormAnalysis/matlab/sean/VabAnalysis/midlineAnalysis/figures/%d/%d.pdf', mvmnt_rep(animal), animal));
end
close(figure);
textprogressbar('done');

%% intensity registration
I = bin410(:,:,1);

raw.i410 = dlmread(fullfile(DIRECTORY, INT_410_FNAME), '', 1, 1);
raw.i470 = dlmread(fullfile(DIRECTORY, INT_470_FNAME), '', 1, 1);
raw.sq410 = ssquare(clip_sj(raw.i410, 1000));
raw.sq470 = ssquare(clip_sj(raw.i470, 1000));
raw.R = raw.sq410 ./ raw.sq470;
raw.OxD = ja_oxd(raw.R);
raw.E = ja_E(raw.OxD);

[reg.fd410, reg.fd470, warp, regInts] =  ...
                ChannelRegister(raw.sq410, raw.sq470, 1000);
reg.i410 = regInts.m410;
reg.i470 = regInts.m470;

reg.R = reg.i410 ./ reg.i470;
reg.OxD = ja_oxd(reg.R);
reg.E = ja_E(reg.OxD);

%% intensity plot
figure;
subplot(2,2,1);
plot(reg.i410(:,1:60));
title('i410 [Animal 1]');
subplot(2,2,2);
plot(reg.i410(:,61:end));
title('i410 [Animal 1]');
subplot(2,2,3);
plot(reg.E(:,1:60));
title('E [Animal 1]');
subplot(2,2,4);
plot(reg.E(:,61:end));
title('E [Animal 2]');

export_fig /Users/sean/Desktop/wormAnalysis/matlab/sean/VabAnalysis/midlineAnalysis/figures/registeredData.pdf

%% intensity error analysis
meanRegi410_1 = mean(reg.i410(:,1:60).').';
meanRegi410_2 = mean(reg.i410(:,61:end).').';

errFromMeani410_1 = reg.i410(:, 1:60) - repmat(meanRegi410_1, [1, 60]);
errFromMeani410_2 = reg.i410(:, 61:end) - repmat(meanRegi410_2, [1, 61]);


meanE_1 = mean(reg.E(:,1:60).').';
meanE_2 = mean(reg.E(:,61:end).').';
errFromMeanE_1 = reg.E(:, 1:60) - repmat(meanE_1, [1, 60]);
errFromMeanE_2 = reg.E(:, 61:end) - repmat(meanE_2, [1, 61]);

% unclipped
figure;
subplot(2,2,1);
plot(errFromMeani410_1); hold on;
yyaxis right;
plot(meanRegi410_1, 'Color', 'k');
title('i410 - mean(i410) [Animal 1]');
subplot(2,2,2);
plot(errFromMeani410_2); hold on;
yyaxis right;
plot(meanRegi410_2, 'Color', 'k');
title('i410 - mean(i410) [Animal 2]');
subplot(2,2,3);
plot(errFromMeanE_1); hold on;
yyaxis right;
plot(meanRegi410_1, 'Color', 'k');
title('subplot 3');
title('E - mean(E) [Animal 1]');
subplot(2,2,4);
plot(errFromMeanE_2); hold on;
yyaxis right;
plot(meanRegi410_2, 'Color', 'k');
title('E - mean(E) [Animal 2]');
export_fig /Users/sean/Desktop/wormAnalysis/matlab/sean/VabAnalysis/midlineAnalysis/figures/errorsFromMeans.pdf


%% clipped E
figure;
subplot(2,2,1);
plot(errFromMeani410_1); hold on;
yyaxis right;
plot(meanRegi410_1, 'Color', 'k');
title('i410 - mean(i410) [Animal 1]');
subplot(2,2,2);
plot(errFromMeani410_2); hold on;
yyaxis right;
plot(meanRegi410_2, 'Color', 'k');
title('i410 - mean(i410) [Animal 2]');
subplot(2,2,3);
plot(errFromMeanE_1);
ylim([-10 10]);
hold on;
yyaxis right;
plot(meanRegi410_1, 'Color', 'k');
title('subplot 3');
title('E - mean(E) [Animal 1]');
subplot(2,2,4);

plot(errFromMeanE_2);
ylim([-10 10]);
hold on;
yyaxis right;
plot(meanRegi410_2, 'Color', 'k');
title('E - mean(E) [Animal 2]');
export_fig /Users/sean/Desktop/wormAnalysis/matlab/sean/VabAnalysis/midlineAnalysis/figures/errorsFromMeans_clippedE.pdf


%% Centroid Distance Curves
I = bwperim(bin410(:,:,1));
stats = regionprops(I);
centroid = stats.Centroid;
[r,c] = size(I);
x = 1:c;
y = 1:r;
[X_,Y_] = meshgrid(x,y) ;
idx = I>0;
figure;
subplot(2,1,1);
scatter(X_(idx), Y_(idx)); hold on;
scatter(centroid(1), centroid(2));
dists = zeros(size(X_(idx), 1));
cmap = colormap;
for i=1:size(X_(idx), 1)
    xs = X_(idx);
    ys = Y_(idx);
    linX = [centroid(1) xs(i)];
    linY = [centroid(2) ys(i)];
    dists(i) = sqrt((xs(i)-centroid(1))^2 + (ys(i) - centroid(2))^2);
    
    line(linX, linY);
end
subplot(2,1,2);
plot(dists(:,1));
export_fig /Users/sean/Desktop/wormAnalysis/matlab/sean/VabAnalysis/midlineAnalysis/figures/centroidToBoundaryLines.pdf

%% Shape spline
x = X_(idx);
y = Y_(idx);

cx = mean(X_(idx));
cy = mean(Y_(idx));
a = atan2(y - cy, x - cx);
[~, order] = sort(a);
x = x(order);
y = y(order);

hold on;
xy = [x y].';
figure;fnplt(cscvn(xy),'r',2);
hold off;
