DIRECTORY = '/Users/sean/Desktop/vab1_2do_05_25_mvmt';
%% Load / Transform data
% mvmnt_rep = readtable(fullfile(DIRECTORY, 'movement.csv'), 'Format', '%d%d%*s');
% mvmnt_rep = table2array(mvmnt_rep(:,2));
COORDS_410_FNAME = 'PA-vab1_2do_05_25_mvmt_410_coords.dat';
COORDS_470_FNAME = 'PA-vab1_2do_05_25_mvmt_470_coords.dat';

INT_410_FNAME = 'PA-vab1_2do_05_25_mvmt_410_intensities.txt';
INT_470_FNAME = 'PA-vab1_2do_05_25_mvmt_470_intensities.txt';

coords410 = loadCoordinates(fullfile(DIRECTORY, COORDS_410_FNAME));
coords470 = loadCoordinates(fullfile(DIRECTORY, COORDS_470_FNAME));

raw.i410 = dlmread(fullfile(DIRECTORY, INT_410_FNAME), '', 1, 1);
raw.i470 = dlmread(fullfile(DIRECTORY, INT_470_FNAME), '', 1, 1);
raw.sq410 = ssquare(clip_sj(raw.i410, 1000));
raw.sq470 = ssquare(clip_sj(raw.i470, 1000));
raw.R = raw.sq410 ./ raw.sq470;
raw.OxD = ja_oxd(raw.R);
raw.E = ja_E(raw.OxD);
%%
[reg.fd410, reg.fd470, warp, regInts] =  ...
                ChannelRegister(raw.sq410, raw.sq470, 1000);
reg.i410 = regInts.m410;
reg.i470 = regInts.m470;

reg.R = reg.i410 ./ reg.i470;
reg.OxD = ja_oxd(reg.R);
reg.E = ja_E(reg.OxD);

%% Plots
mkdir(fullfile(DIRECTORY, 'figures'));
mkdir(fullfile(DIRECTORY, 'figures', 'registration'));
for i = 1:size(reg.i410, 2)
    plot(reg.i410(:,i)); hold on; plot(reg.i470(:,i)); hold off; 
    export_fig(fullfile(DIRECTORY, 'figures', 'registration', strcat(num2str(i), '.pdf')))
end

%% Warp Function Error
line = linspace(1,100,1000);
for i = 1:size(reg.i410, 2)
    err(:,i) = abs((warp(:,i).' - line).');
    sum_err(i) = sum(err(:,i));
end

%% mean E values
% mean(reg.E).';

%% Midline Analysis
coords_410_xy = table2array(PAvab12do0525mvmt410coords);
coords_470_xy = table2array(PAvab12do0525mvmt470coords);

coords410.x = coords_410_xy(:,1:2:end);
coords410.y = coords_410_xy(:,2:2:end);

coords470.x = coords_470_xy(:,1:2:end);
coords470.y = coords_470_xy(:,2:2:end);