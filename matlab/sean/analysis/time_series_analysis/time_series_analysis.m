imageFile = "/Users/sean/code/wormAnalysis/data/cata_time_series/120827_HD233_WT_2doF_5mMTB_ALLplts/120827_HD233_WT_2do_plt1/0_Stacks 120827_HD233_WT_2do_plt1/120827_HD233_WT_2doF_PreLev_plt1_01.tif";
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

% Flip AP
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

%% Focus Analysis
bbox_for_rot = [134 105 80 50]; % x0 y0 w h
focus410 = zeros(nAnimals, 1);
focus470 = zeros(nAnimals, 1);
for i=1:nAnimals
    focus410(i) = fmeasure(im410(:,:,i), 'CURV', bbox_for_rot);
    focus470(i) = fmeasure(im470(:,:,i), 'CURV', bbox_for_rot);
end
scatter(mean_i410, blurAnnotations.Value);
% scatter(focus410,blurAnnotations.Value);

%%
indexer_table = readtable("/Users/sean/code/wormAnalysis/data/cata_time_series/120827_HD233_WT_2doF_5mMTB_ALLplts/120827_HD233_WT_2do_plt1/indexer.csv");

%%
strain_frame_mapper = table(...
    'Size', [1936 2], 'VariableTypes', {'string', 'uint16'},...
    'VariableNames', {'Strain', 'Frame'});
for i=1:height(indexer_table)
    strain = indexer_table.Strain{i};
    start_ = indexer_table.Start_Animal(i);
    end_ = indexer_table.End_Animal(i);
    
    strain_frame_mapper.Frame(start_:end_) = start_:end_;
    strain_frame_mapper.Strain(start_:end_) = strain;
end

writetable(strain_frame_mapper, '~/Desktop/strain_mapper.csv');

%%