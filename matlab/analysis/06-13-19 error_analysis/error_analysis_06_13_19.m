splitImages = splitImageStack(loadTiffImageStack(fullfile("data", "2017_02_22-HD233_SAY47.tif")), "TL/470/410/470/410");

%
% [i410_1_cata, i470_1_cata] = pipelineCata(splitImages.imTL_1, splitImages.im410_1, splitImages.im470_1);
% [i410_1_new, i470_1_new] = pipelineTwoMidlinesTwoMasks(splitImages.imTL_1, splitImages.im410_1, splitImages.im470_1);
% 
% [i410_2_cata, i470_2_cata] = pipelineCata(splitImages.imTL_1, splitImages.im410_2, splitImages.im470_2);
% [i410_2_new, i470_2_new] = pipelineTwoMidlinesTwoMasks(splitImages.imTL_1, splitImages.im410_2, splitImages.im470_2);
% 
% R_1_cata = i410_1_cata ./ i410_2_cata;
% R_1_new = i410_1_new ./ i410_2_new;

%
% Generate a new "synthetic" dataset using only i410 data
[centroids, orientations] = calcCentroidsAndOrientations(segmentPharynx(splitImages.im410_1,0,2000));
orientations = deg2rad(-orientations);
%%
minShift = 0;
maxShift = 3;
nShifts = 5;

unitXY = [cos(orientations) sin(orientations)]; % [x y]
shiftMagnitudes = linspace(minShift, maxShift, nShifts);
unitXYZ = bsxfun(@times,unitXY,reshape(shiftMagnitudes,1,1,[])); % [animal (x1 y2) shiftIdx]

unshiftedIm410 = zeros(size(splitImages.im410_1, 1), size(splitImages.im410_1, 2), size(splitImages.im410_1, 3) * nShifts);
shiftedIm410 = unshiftedIm410;

shiftMagnitudesRep = repmat(shiftMagnitudes, [1 size(splitImages.im410_1, 3)]);
for i=1:size(shiftedIm410, 3)
    animalIndex = ceil(i/nShifts);
    zIndex = mod(i-1, nShifts)+1;
    unshiftedIm410(:,:,i) = splitImages.im410_1(:,:,animalIndex);
    shiftedIm410(:,:,i) = imtranslate(splitImages.im410_1(:,:,animalIndex), [unitXYZ(animalIndex, 1, zIndex), unitXYZ(animalIndex, 2, zIndex)]);
end

imE = double(unshiftedIm410) - double(shiftedIm410);
% implay(imE);
options.overwrite = true;
saveastiff(single(imE), '~/Desktop/imE_bw.tiff', options);

[shifted_i410_1_cata, shifted_i410_2_cata] = pipelineCata(splitImages.imTL_1, unshiftedIm410, shiftedIm410);
[shifted_i410_1_new, shifted_i410_2_new] = pipelineTwoMidlinesTwoMasks(splitImages.imTL_1, unshiftedIm410, shiftedIm410);
% 
shifted_percent_e_cata = abs(shifted_i410_1_cata - shifted_i410_2_cata) ./ ((shifted_i410_1_cata + shifted_i410_2_cata) / 2);
shifted_percent_e_new = abs(shifted_i410_1_new - shifted_i410_2_new) ./ ((shifted_i410_1_new + shifted_i410_2_new) / 2);

% % plot
% 
% mean_shifted_e_cata = mean(shifted_percent_e_cata, 1);
% mean_shifted_e_new = mean(shifted_percent_e_new, 1);
% 
% plotMultiplePharynxData({shifted_percent_e_cata, shifted_percent_e_new}, {'cata', 'new'});
% 
% figure;
% scatter(abs(shiftMagnitudesRep), mean_shifted_e_cata);
% title('Cata');
% figure;
% scatter(abs(shiftMagnitudesRep), mean_shifted_e_new);
% title('New');
% 
% %
% t = table(mean_shifted_e_cata.', mean_shifted_e_new.', shiftMagnitudesRep.');
% writetable(t, '~/Desktop/shifted_percent_errors.csv');

%% Contours
contour_new = {};
contour_cata = {};
for i=1:nShifts
    contour_new{end+1} = 100*shifted_percent_e_new(:,i:nShifts:end);
    contour_cata{end+1} = 100*shifted_percent_e_cata(:,i:nShifts:end);
end

figure;
ax = subplot(2,1,1);
plotMultiplePharynxData(contour_cata, sprintfc('%0.1f',shiftMagnitudes), [0 28] , ax);
title(ax, 'Old Pipeline');

ax = subplot(2,1,2);
plotMultiplePharynxData(contour_new, sprintfc('%0.1f',shiftMagnitudes), [0 7] ,ax);
title(ax, 'New Pipeline');