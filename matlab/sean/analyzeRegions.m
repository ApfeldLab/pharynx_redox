% Load an experiment
% You should only change e_name per experiment

root_exp_dir  = '/Volumes/dept-shares/Image Analysis/data/';
e_name = '2018_06_07_SAY98_HD233';

%%
% Channel registration
e = Experiment(fullfile(root_exp_dir, e_name));
e.registerChannels();
save(fullfile(root_exp_dir, e_name, e_name), 'e');

% Save region data
writetable([e.metadata e.reg.regions.all], fullfile(root_exp_dir, e_name, strcat(e_name, '-alldatawithregeions.csv')));

%% Plots
