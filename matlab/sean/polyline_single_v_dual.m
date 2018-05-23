% Question:
% On the images labelled "0", how similar is using two polylines (one per
% channel) as compared to using a single polyline (taken on 410) for both
% channels?

%% Load Data

% im_410 = tiffread2('../data/proc_data/all_410_PA.tif');
% im_470 = tiffread2('../data/proc_data/all_470_PA.tif');

mvmt_labels = csvread('../data/raw_data_1/labels.csv');

p_410_m_410 = csvread('../experiments/single_v_dual_poly/poly_410_measure_410_sean.csv', 1, 1);
p_410_m_470 = csvread('../experiments/single_v_dual_poly/poly_410_measure_470_sean.csv', 1, 1);

p_470_m_410 = csvread('../experiments/single_v_dual_poly/poly_470_measure_410_sean.csv', 1, 1);
p_470_m_470 = csvread('../experiments/single_v_dual_poly/poly_470_measure_470_sean.csv', 1, 1);

%% Plot
plot(p_410_m_470)