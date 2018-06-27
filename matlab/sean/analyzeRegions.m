% Load an experiment
% You should only change e_name per experiment

root_exp_dir  = '/Users/sean/Desktop/';
exp_name = '2018_06_14_SAY98_HD233';

%%

experiment_directory = fullfile(root_exp_dir, exp_name);

% Channel registration
e = Experiment(experiment_directory);
e.registerChannels();
save(fullfile(experiment_directory, exp_name), 'e');

% Save region data
writetable(e.reg.regions.all, fullfile(root_exp_dir, exp_name, strcat(exp_name, '_region_data.csv')));

%% Plots
