% Load metadata right after ImageJ Analysis

% Get experiment name, defaults to directory name, ultimately asks user

directory = '/Users/sean/Desktop/vab1_2do_05_25_mvmt/';

EXP_NAME = 'vab1_2do_05_25_mvmt';

% Load ImageJ Analysis data
coords_410 = loadDataFile(directory, strcat('PA-', EXP_NAME ,'_410_coords.txt'));
coords_470 = loadDataFile(directory, strcat('PA-', EXP_NAME ,'_470_coords.txt'));
intensity_410 = loadDataFile(directory, strcat('PA-', EXP_NAME ,'_410_intensities.txt'));
intensity_470 = loadDataFile(directory, strcat('PA-', EXP_NAME ,'_470_intensities.txt'));

% Trim tables, because ImageJ outputs row numbers etc.
coords_410 = coords_410(:, 2:end);
coords_470 = coords_470(:, 2:end);
intensity_410 = table2array(intensity_410(2:end, 2:end));
intensity_470 = table2array(intensity_470(2:end, 2:end));

% Length-normalize to 100xn vectors
% where n = number of animals
sq_intensity_410 = ssquare(intensity_410);
sq_intensity_470 = ssquare(intensity_470);

[reg410_FD, reg470_FD, warp470, resampled_intensity, fdObjs] = ChannelRegister(sq_intensity_410, sq_intensity_470, 1000);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkIfDir(directory)
    if ~isfolder(directory)
        error('Error. The supplied directory is not a directory');
    end
end

function dataTable = loadDataFile(directory, globPattern)
    checkIfDir(directory);
    
    data_filename = dir(fullfile(directory, globPattern));
    if size(data_filename, 1) > 1
        error('Error. More than 1 measurement file in given directory. Delete one and try again');
    else
        data_filename = fullfile(data_filename(1).folder, data_filename(1).name);
    end
    
    dataTable = readtable(data_filename);
end