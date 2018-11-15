rawImgFilePath = "/Users/sean/Documents/worm stuff/vab1_2do_05_25_mvmt_new/vab1_2do_05_25_mvmt.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);
allImgs = loadTiffImageStack(rawImgFilePath);

imTL = allImgs(:, :, 1:3:end);
im410_raw = allImgs(:, :, 2:3:end);
im470_raw = allImgs(:, :, 3:3:end);
nAnimals = size(im410_raw, 3);

im_r_410_470 = im410_raw ./ im470_raw;

mvmnt = csvread(fullfile(rawImgDir, 'movement.csv'));

im410 = subtractMedians(im410_raw);
im470 = subtractMedians(im470_raw);

seg410 = segmentPharynx(im410, 0);
seg470 = segmentPharynx(im470, 0);

segLRBounds = double(getLeftRightBounds(seg410));
%%
midlines410 = calculateMidlines(imTL, seg410);
midlines470 = calculateMidlines(imTL, seg470);

%%
INTERP_METHOD = 'bilinear';
PROFILE_LENGTH = 1000;
imWidth = size(im410, 2);

i410_raw = measureIntensityAlongMidlines(im410, midlines410, segLRBounds, PROFILE_LENGTH, INTERP_METHOD);
i470_raw = measureIntensityAlongMidlines(im470, midlines470, segLRBounds, PROFILE_LENGTH, INTERP_METHOD);

[trimmed410, unscaledBounds410, scaledBounds410] = trimProfile(i410_raw, segLRBounds);
[trimmed470, unscaledBounds470, scaledBounds470] = trimProfile(i470_raw, segLRBounds);

i410 = ssquare(trimmed410);
i470 = ssquare(trimmed470);

%% Functionalize
% fd_410 = makeIntensityFdObj(i410);
% fd_470 = makeIntensityFdObj(i470);

[reg410_FD, reg470_FD, warp470, resampled_intensity, fdObjs] = ...
    ChannelRegister(i410, i470, 100);

%% Midline Stuff
animal_number = 10;
xs410 = linspace(scaledBounds410(animal_number, 1), scaledBounds410(animal_number, 2), 100);
ys410 = feval(midlines410{animal_number}, xs410);
ys470 = feval(midlines470{animal_number}, xs410);

figure;
imshow(im410_raw(:,:,animal_number), []); hold on;
% plot(midlines410{animal_number});
plot(xs410, ys410); hold on;
plot(xs410, ys470);

hold off;

%% Midline Warping

animal_number = 75;
PROF_SAMPLE_LEN = 5000;

xs410 = linspace(scaledBounds410(animal_number, 1), scaledBounds410(animal_number, 2), PROF_SAMPLE_LEN);
xs470 = linspace(scaledBounds470(animal_number, 1), scaledBounds470(animal_number, 2), PROF_SAMPLE_LEN);

ys410 = feval(midlines410{animal_number}, xs410);
ys470 = feval(midlines470{animal_number}, xs470);

figure;
plot(xs410, ys410); hold on;
plot(xs470, ys470); hold off;

s = linspace(1, 100, PROF_SAMPLE_LEN);
warp_s = eval_fd(s, fdObjs(animal_number).warpFD);
warp_s = warp_s(:,2);

warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);

warped_ys470 = ys470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));
warped_xs470 = xs470(min(round(warp_s_scaled), PROF_SAMPLE_LEN));

w = eval_fd(s, fdObjs(animal_number).warpFD);
w = w(:,2);

% Midline warp plot
figure;
% subplot(2, 1, 1);
plot(xs410, ys410); hold on;
plot(xs470, ys470);
plot(warped_xs470, warped_ys470, 'LineWidth', 1); 
yyaxis right;
plot(s,s,'g--');
plot(s, w, 'g');
legend('410', '470', 'w470'); hold off;

% distance plot
figure;
Axy = [xs410; ys410.'];
Bxy = [warped_xs470; warped_ys470.'];

d = sqrt(sum((Axy-Bxy).^2,1));
plot(s, d);
ylim([0 3]);
yyaxis right;
plot(s,s,'r--'); hold on;
plot(s, w, 'r');hold off;

%% Useful warp plot
animal_number = 31;
figure;
% subplot(2, 1, 2);
s = linspace(1, 100, PROF_SAMPLE_LEN);
reg410ys = eval_fd(s, reg410_FD(animal_number)); 
% reg410ys = reg410ys(:,2);
reg470ys = eval_fd(s, reg470_FD(animal_number)); 
% reg470ys = reg470ys(:,2);
orig410ys = eval_fd(s, fdObjs(animal_number).origFD); 
% orig410ys = orig410ys(:,2);

plot(s, reg410ys, 'b-'); hold on;
plot(s, reg470ys, 'r-');
plot(s, orig410ys, 'k--');
yyaxis right;
plot(s,s);
w = eval_fd(s, fdObjs(animal_number).warpFD);
w = w(:,2);
plot(s,s,'g--');
plot(s, w, 'g');
legend('r410', 'r470', 'o470', 'nowarp', 'warp');
hold off;

%% Multiple midline warp distance plots
PROF_SAMPLE_LEN = 5000;
my_rgb = cbrewer('qual', 'Set1', 5);
% lw = [.05 4];
figure; hold on;
for animal_number=1:nAnimals
    xs410 = linspace(scaledBounds410(animal_number, 1), scaledBounds410(animal_number, 2), PROF_SAMPLE_LEN);
    ys410 = feval(midlines410{animal_number}, xs410);
    ys470 = feval(midlines470{animal_number}, xs410);

    s = linspace(1, 100, PROF_SAMPLE_LEN);
    warp_s = eval_fd(s, fdObjs(animal_number).warpFD);
    warp_s = warp_s(:,2);
    warp_s_scaled = rescale(warp_s, 1, PROF_SAMPLE_LEN);

    warped_ys470 = ys470(min(floor(warp_s_scaled), PROF_SAMPLE_LEN));
    warped_xs470 = xs410(min(floor(warp_s_scaled), PROF_SAMPLE_LEN));

    w = eval_fd(s, fdObjs(animal_number).warpFD);
    w = w(:,2);

    % distance plot
    
    Axy = [xs410; ys410.'];
    Bxy = [warped_xs470; warped_ys470.'];

%     d = sqrt(sum((Axy-Bxy).^2,1));
    d = xs410 - warped_xs470;
    plot(s, d, 'Color', my_rgb(mvmnt(animal_number), :));
    ylim([0 3]);
end
hold off;
