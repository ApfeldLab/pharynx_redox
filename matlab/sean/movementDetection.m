%% Load Data
rawImgFilePath = "/Users/sean/Documents/worm stuff/vab1_2do_05_25_mvmt_new/vab1_2do_05_25_mvmt.tif";
[rawImgDir, ~, ~]= fileparts(rawImgFilePath);

allImgs = loadTiffImageStack(rawImgFilePath);

imTL = allImgs(:, :, 1:3:end);
im410 = double(allImgs(:, :, 2:3:end));
im470 = double(allImgs(:, :, 3:3:end));
nAnimals = size(im410, 3);

im410(im410 == 0) = 1;
im470(im470 == 0) = 1;

im_r_410_470 = im410 ./ im470;

seg410 = segmentPharynx(im410, 0);
seg470 = segmentPharynx(im470, 0);

midlines410 = calculateMidlines(imTL, seg410);
midlines470 = calculateMidlines(imTL, seg470);

INTERP_METHOD = 'bilinear';
PROFILE_LENGTH = 1000;
i410 = measureIntensityAlongMidlines(im410, midlines410, PROFILE_LENGTH, INTERP_METHOD);
i470 = measureIntensityAlongMidlines(im470, midlines410, PROFILE_LENGTH, INTERP_METHOD);
i_r_410_470 = measureIntensityAlongMidlines(im_r_410_470, midlines410, PROFILE_LENGTH, INTERP_METHOD);

mvmnt = csvread(fullfile(rawImgDir, 'movement.csv'));
nAnimals = size(im410, 3);
%%
[i410_trimmed, bounds] = trimProfile(i410);
sq_i410 = ssquare(i410_trimmed);

%%
i_r_410_470_trimmed = i_r_410_470;
for i=1:nAnimals
    i_r_410_470_trimmed(1:bounds(i,1),i) = 0;
    i_r_410_470_trimmed(bounds(i,2):end,i) = 0;
end

sq_r_410_470 = ssquare(i_r_410_470_trimmed);
sq_r_410_470_squared = sq_r_410_470 .^ 2;

any_mvmt = mvmnt > 0;

% plotStrainMeans(num2str(any_mvmt), sq_r_410_470, "Ratio Movement Means");
%%
n = 10; % every n values

no_mvmt = sq_r_410_470(:,mvmnt == 0);
some_mvmt = sq_r_410_470(:,mvmnt ~= 0);

means = zeros(n,2); % None, Some
stds = zeros(n,2);  % None, Some
variances = zeros(n,2); % None, Some
for i=1:n
    s = ((i-1) * n) + 1;
    e = i * n;
    
    none = no_mvmt(s:e, :);
    some = some_mvmt(s:e, :);
    
    means(i, :) = [mean(none(:)) mean(some(:))];
    stds(i, :) = [std(none(:)) std(some(:))];
    variances(i, :) = [var(none(:)) var(some(:))];
end

bar(stds)
legend('none', 'some');

%%
unique_movements = unique(mvmnt);
nMovmentCategories = size(unique_movements, 1);
chunkSize = 10;
means_block = zeros(uint8(size(sq_r_410_470, 1) / chunkSize), nMovmentCategories);
stds_block = zeros(uint8(size(sq_r_410_470, 1) / chunkSize), nMovmentCategories);
all_stds = zeros(uint8(size(sq_r_410_470, 1) / chunkSize), 1);
for i=1:nMovmentCategories
    data_ = sq_r_410_470(:,mvmnt == i);
    if ~isempty(data_)
        means_block(:,i) = blockproc(data_, [chunkSize size(data_, 2)], @(bs) mean(bs.data(:)));
        stds_block(:,i) = blockproc(data_, [chunkSize size(data_, 2)], @(bs) std(bs.data(:)));
    end
end
bar(stds_block);
legend(num2str(unique_movements));

%% PCA on StdDev
chunkSize = 10;
all_stds =  blockproc(sq_r_410_470, [chunkSize 1], @(bs) std(bs.data(:))).';
[coeff, score, latent, tsquared, explained] = pca(all_stds, 'NumComponents', 3);

%% 3d PCA scatter
figure; 
scatter3(score(mvmnt==0,1), score(mvmnt==0,2), score(mvmnt==0,3)); hold on;
scatter3(score(mvmnt==1,1), score(mvmnt==1,2), score(mvmnt==1,3)); hold on;
scatter3(score(mvmnt==2,1), score(mvmnt==2,2), score(mvmnt==2,3)); hold on;
scatter3(score(mvmnt==3,1), score(mvmnt==3,2), score(mvmnt==3,3)); hold on;
legend('0', '1', '2', '3');
hold off;

%% 2d PCA scatter
figure; 
scatter(score(mvmnt==0,1), score(mvmnt==0,2)); hold on;
scatter(score(mvmnt==1,1), score(mvmnt==1,2)); hold on;
scatter(score(mvmnt==2,1), score(mvmnt==2,2)); hold on;
scatter(score(mvmnt==3,1), score(mvmnt==3,2)); hold on;
legend('0', '1', '2', '3');
hold off;

%% Classifier stuff
stds_and_labels = horzcat([all_stds mvmnt ~= 0]);

%%
i = 4;
xs = 1:size(im410, 3);
imshow(im_r_410_470(:,:,i), []); hold on;
plot(xs, feval(midlines410{i}, xs), 'color', 'r', 'LineWidth', 2); hold on;
plot(xs, feval(midlines470{i}, xs), 'color', 'b', 'LineWidth', 2); hold off;
legend('410', '470');