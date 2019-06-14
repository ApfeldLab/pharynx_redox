% Load Data
rawImgFilePath = "/Users/sean/Desktop/midlineMovment/WT high mvmt 11-19-18.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
allImgs = loadTiffImageStack(rawImgFilePath);

% Frame 1 -- animal #1, transmitted light
% Frames 2-127 -- animal #1, 470/410
% Frames 128-241 -- animal #1, 410/410
% Frame 242 -- animal #2, transmitted light
% Frames 243-364 -- animal #2, 470/410
% Frames 365-488 -- animal #2, 410/410
%%

im470_raw = allImgs(:, :, 128:2:241);
im470_raw = cat(3, im470_raw, allImgs(:, :, 365:2:488));
im410_raw = allImgs(:, :, 129:2:241);
im410_raw = cat(3, im410_raw, allImgs(:, :, 366:2:488));

nAnimals = size(im410_raw, 3);
imTL = repmat(allImgs(:, :, 1), int8([1,1,numel(128:2:241)]));
imTL = cat(3, imTL, repmat(allImgs(:, :, 242), int8([1,1,numel(365:2:488)])));
im_r_410_470 = im410_raw ./ im470_raw;
mvmnt = csvread('/Users/sean/Desktop/midlineMovment/movement.csv');

%% Process Data
im410 = subtractMedians(im410_raw);
im470 = subtractMedians(im470_raw);

seg410 = segmentPharynx(im410, 0);
seg470 = segmentPharynx(im470, 0);

segLRBounds = double(getLeftRightBounds(seg410));

midlines410 = calculateMidlines(imTL, seg410);
midlines470 = calculateMidlines(imTL, seg470);

INTERP_METHOD = 'bilinear';
PROFILE_LENGTH = 1000;
imWidth = size(im410, 2);

i410_raw = measureIntensityAlongMidlines(im410, midlines410, segLRBounds, PROFILE_LENGTH, INTERP_METHOD);
i470_raw = measureIntensityAlongMidlines(im470, midlines470, segLRBounds, PROFILE_LENGTH, INTERP_METHOD);

[trimmed410, unscaledBounds410, scaledBounds410] = trimProfile(i410_raw, segLRBounds);
[trimmed470, unscaledBounds470, scaledBounds470] = trimProfile(i470_raw, segLRBounds);

i410 = ssquare(trimmed410);
i470 = ssquare(trimmed470);

[reg410_FD, reg470_FD, warp470, resampled_intensity, fdObjs] = ...
    ChannelRegister(i410, i470, 100);

%% Main Plot
for i=1:nAnimals
    PROF_SAMPLE_LEN = 100;

    xs410 = linspace(scaledBounds410(i, 1), scaledBounds410(i, 2), PROF_SAMPLE_LEN);
    xs470 = linspace(scaledBounds470(i, 1), scaledBounds470(i, 2), PROF_SAMPLE_LEN);

    ys410 = feval(midlines410{i}, xs410);
    ys470 = feval(midlines470{i}, xs470);

    s = linspace(1, 100, PROF_SAMPLE_LEN);
    warp_s = eval_fd(s, fdObjs(i).warpFD);
    warp_s = warp_s(:,2);

    warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);

    warped_ys470 = ys470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
    warped_xs470 = xs470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));

    w = eval_fd(s, fdObjs(i).warpFD);
    w = w(:,2);
    Axy = [xs410; ys410.'];
    Bxy = [warped_xs470; warped_ys470.'];

    d = sqrt(sum((Axy-Bxy).^2,1));
    subplot(2,2,1);
    plot(s, d);
    ylim([0 2]);
    
    subplot(2,2,2);
    bbox = regionprops(seg410(:,:,i), 'BoundingBox');
    bbox = bbox.BoundingBox;
    xLeft = bbox(1);
    xRight = xLeft + bbox(3);
    yTop = bbox(2);
    yBot = yTop + bbox(4);

    imagesc(im_r_410_470(:,:,i)); hold on;
    % plot(xs410, ys410);
    % plot(warped_xs470, warped_ys470);
    rate = 50;
    region_names = fieldnames(Constants.regions);
    for j=1:numel(region_names) - 2
        region_name = region_names{j};
        region = Constants.regions.(region_name);

        xs410 = linspace(scaledBounds410(i, 1), scaledBounds410(i, 2), PROF_SAMPLE_LEN);
        xs470 = linspace(scaledBounds470(i, 1), scaledBounds470(i, 2), PROF_SAMPLE_LEN);

        ys410 = feval(midlines410{i}, xs410);
        ys470 = feval(midlines470{i}, xs470);

        s = linspace(1, 100, PROF_SAMPLE_LEN);
        warp_s = eval_fd(s, fdObjs(i).warpFD);
        warp_s = warp_s(:,2);

        warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);

        warped_ys470 = ys470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
        warped_xs470 = xs470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));

        plot(warped_xs470(1, region(1):region(2)), warped_ys470(region(1):region(2), 1), 'LineWidth', 5);
%     legend(region_names(1:end-2));
    end
    
    plot([downsample(xs410, rate); downsample(warped_xs470, rate)], [downsample(ys410.', rate); downsample(warped_ys470.', rate)], 'Color', 'red');
    xlim([xLeft xRight]);
    ylim([yTop yBot]);
    hold off;


    subplot(2,2,[3 4]);
    set(gcf, 'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
    xs = linspace(1,100,200);
    plot(xs, eval_fd(xs,reg410_FD(i)), 'Color', 'blue'); hold on;
    plot(xs, eval_fd(xs,reg470_FD(i)), 'Color', 'red');
    ylim([0 15000]);
    yyaxis right;
    plot(1:100,1:100, 'Color', 'black');
    plot(warp470(:,i), 'Color', 'red');hold off;
    xlim([1 100]);
    ylim([1 100]);
    
%     set(gcf, 'Position', get(0, 'Screensize'));
%     export_fig(sprintf('/Users/sean/Desktop/midlineMovment/%d', i), '-q101', '-a1', '-pdf');
%     clf;
end

%% Region Ratio Table
t = regionMeans(resampled_intensity.m410 ./ resampled_intensity.m470);
t.Movement = mvmnt;
t.Animal = cat(2, repmat(1, [1 numel(128:2:241)]), repmat(2, [1 numel(365:2:488)])).';
writetable(t, '/Users/sean/Desktop/midlineMovment/regionMeans.csv');

%% Unregistered Region Ratio Table
t_unreg = regionMeans(i410 ./ i470);
t_unreg.Movement = mvmnt;
t_unreg.Animal = cat(2, repmat(1, [1 numel(128:2:241)]), repmat(2, [1 numel(365:2:488)])).';
writetable(t_unreg, '/Users/sean/Desktop/midlineMovment/regionMeans_unreg.csv');

%% Save images
imwrite(im_r_410_470(:,:,1),'/Users/sean/Desktop/midlineMovment/r410_410.tif');
% for i=2:nAnimals
%     imwrite(im_r_410_470(:,:,i),'/Users/sean/Desktop/midlineMovment/r410_410.tif', 'WriteMode', 'append');
% end

%% Create Dynamic Time Warp Distance Matrix
PROF_SAMPLE_LEN = 100;

dtw_d = zeros(nAnimals, PROF_SAMPLE_LEN);
dtw_dx = zeros(nAnimals, PROF_SAMPLE_LEN);
dtw_dy = zeros(nAnimals, PROF_SAMPLE_LEN);

% for i=1:nAnimals
%     xs410 = linspace(scaledBounds410(i, 1), scaledBounds410(i, 2), PROF_SAMPLE_LEN);
%     xs470 = linspace(scaledBounds470(i, 1), scaledBounds470(i, 2), PROF_SAMPLE_LEN);
% 
%     ys410 = feval(midlines410{i}, xs410);
%     ys470 = feval(midlines470{i}, xs470);
% 
%     s = linspace(1, 100, PROF_SAMPLE_LEN);
%     warp_s = eval_fd(s, fdObjs(i).warpFD);
%     warp_s = warp_s(:,2);
% 
%     warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);
% 
%     warped_ys470 = ys470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
%     warped_xs470 = xs470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
%     
%     
% 
%     w = eval_fd(s, fdObjs(i).warpFD);
%     w = w(:,2);
%     Axy = [xs410; ys410.'];
%     Bxy = [warped_xs470; warped_ys470.'];
% 
%     dtw_d(i, :) = sqrt(sum((Axy-Bxy).^2,1));
%     dtw_dx(i, :) = xs410 - warped_xs470;
%     dtw_dy(i, :) = ys410 - warped_ys470;
% end

%% Tables for the DTW distances
t_dtw_d = regionMeans(dtw_d);
t_dtw_dx = regionMeans(dtw_dx);
t_dtw_dy = regionMeans(dtw_dy);

writetable(t_dtw_d, '/Users/sean/Desktop/midlineMovment/regionMeans_dtw_d.csv');
writetable(t_dtw_dx, '/Users/sean/Desktop/midlineMovment/regionMeans_dtw_dx.csv');
writetable(t_dtw_dy, '/Users/sean/Desktop/midlineMovment/regionMeans_dtw_dy.csv');

%% Other plot
% for i=1:nAnimals
% imagesc(im410(:,:,i)); hold on;
% region_names = fieldnames(Constants.regions);
% for j=1:numel(region_names) - 2
%     region_name = region_names{j};
%     region = Constants.regions.(region_name);
%     
%     xs410 = linspace(scaledBounds410(i, 1), scaledBounds410(i, 2), PROF_SAMPLE_LEN);
%     xs470 = linspace(scaledBounds470(i, 1), scaledBounds470(i, 2), PROF_SAMPLE_LEN);
% 
%     ys410 = feval(midlines410{i}, xs410);
%     ys470 = feval(midlines470{i}, xs470);
% 
%     s = linspace(1, 100, PROF_SAMPLE_LEN);
%     warp_s = eval_fd(s, fdObjs(i).warpFD);
%     warp_s = warp_s(:,2);
% 
%     warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);
% 
%     warped_ys470 = ys470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
%     warped_xs470 = xs470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
%    
%     plot(warped_xs470(1, region(1):region(2)), warped_ys470(region(1):region(2), 1), 'LineWidth', 5);
%     legend(region_names(1:end-2));
% end
% hold off;
% set(gcf, 'Position', get(0, 'Screensize'));
% export_fig(sprintf('/Users/sean/Desktop/midlineMovment/regions/%d', i), '-q101', '-a1', '-pdf');
% clf;
% end

%%
csvwrite('/Users/sean/Desktop/midlineMovment/i410.csv', i410.')
csvwrite('/Users/sean/Desktop/midlineMovment/i470.csv', i470.')

%% intensity Covariance plot
figure;
a = reg410_FD(1:50);
b = reg470_FD(1);

var_bifd = var_fd(a, a);
s = 1:100;
var_mat = eval_bifd(s, s, var_bifd);
surf(s, s, var_mat);