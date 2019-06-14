% Blur Analysis

% 
%% LOAD DATA

% if .mat file exists 
%   load .mat file
% else
%   do all this stuff V
%   save to .mat

imageFile = "/Users/sean/code/wormAnalysis/data/cata_time_series/120827_HD233_WT_2doF_5mMTB_ALLplts/120827_HD233_WT_2do_plt1/0_Stacks 120827_HD233_WT_2do_plt1/120827_HD233_WT_2doF_5mM_plt1_01.tif";

allImages = loadTiffImageStack(imageFile);

im410 = allImages(:,:,1:2:end);
im470 = allImages(:,:,2:2:end);

im410 = subtractMedians(im410);
im470 = subtractMedians(im470);

nAnimals = size(im410, 3);

t = 800;

seg410 = im410;
seg410(seg410 <  t) = 0;
seg410(seg410 >= t) = 1;

seg470 = im470;
seg470(seg470 <  t) = 0;
seg470(seg470 >= t) = 1;
%
[rot_im410, rot_seg410] = rotatePharynx(im410, seg410);
[rot_im470, rot_seg470] = rotatePharynx(im470, seg470);
%
midlines410 = calculateMidlinesNoTL(rot_seg410);
midlines470 = calculateMidlinesNoTL(rot_seg470);

N_DATA_POINTS = 1000;
%
[i410, i410_raw] = measureAndTrim(rot_im410, midlines410, 500, N_DATA_POINTS);
[i470, i470_raw] = measureAndTrim(rot_im470, midlines470, 500, N_DATA_POINTS);

unflipped_i410 = i410;
unflipped_i470 = i470;


% Flip Anterior-Posterior
shouldFlipAP = calcShouldFlipAP(i410);
for i=1:size(i410,2)
    if shouldFlipAP(i)
        % Data
        i410(:,i) = flip(i410(:,i));
        i470(:,i) = flip(i470(:,i));
        i410_raw(:,i) = flip(i410_raw(:,i));
        i470_raw(:,i) = flip(i470_raw(:,i));
        % Images
        im410(:,:,i) = fliplr(im410(:,:,1));
        im470(:,:,i) = fliplr(im470(:,:,1));
        
    end
end
plot(i410_raw);

%% BLUR ANALYSIS

% Look up:
%   confusion matrix
%   AUC

%% Helper functions
