% A set of helper functions to load and coalesce experiment metadata from 
% the experiment directory.
%
% ASSUMPTIONS:
%   - You have already run ImageJ analysis.
%
% REQUIRED FILES:
%   - EXP_NAME_measurements.csv
%       * This contains image-level and whole-pharynx level intensity and
%       morphological measurements.
%   - EXP_NAME_mvmt.csv
%       * This contains a single column with the movement calls on the
%       ratio images.
%   - EXP_NAME_410-PA_coords.txt
%       * 410nm polyline coordinates
%   - EXP_NAME_470-PA_coords.txt
%       * 470nm polyline coordinates
%   - EXP_NAME_410-PA_intensities.txt
%       * 410nm polyline intensity measurements
%   - EXP_NAME_470-PA_intensities.txt
%       * 470nm polyline intensity measurements
%
% TODO: suppress warnings
function metadata = loadMetadataFromDirectory(directory)
    experimentID = getExperimentID(directory);
%     measurementTable = loadMeasurements(directory);
%     movementRep = loadMvmtRep(directory);
%     metadata = coalesce(experimentID, measurementTable, movementRep);
    metadata = experimentID;
end

function metadata = coalesce(experimentID, measTable, mvmtTable)
    metadata = horzcat(measTable, mvmtTable);
    expIDs = cell(size(metadata,1), 1);
    expIDs(:) = {experimentID};
    metadata.ExperimentID = expIDs;
end

function experimentID = getExperimentID(directory)
    checkIfDir(directory);
    
    [~,EXP_NAME,~] = fileparts(directory);
    experimentID = input(sprintf('Enter an Experiment ID [%s]:', EXP_NAME), 's');
    if strlength(experimentID) == 0
        experimentID = EXP_NAME;
    end
end

function measurementTable = loadMeasurements(directory)
    measurementTable = loadDataFile(directory, '*_measurements.csv');
end

function mvmtRepDataTable = loadMvmtRep(directory)
    mvmtRepDataTable = loadDataFile(directory, '*_mvmt.csv');
    mvmtRepDataTable.Properties.VariableNames = {'Mvmt_Rep'};
end

function dataTable = loadDataFile(directory, globPattern)
    checkIfDir(directory);
    
    data_filename = dir(fullfile(directory, globPattern));
    if size(data_filename, 1) > 1
        error('Error. More than 1 measurement file in given directory. Delete one and try again');
    else
        data_filename = fullfile(data_filename(1).folder, data_filename(1).name);
    end
    
    dataTable = readtable(data_filename, 'Delimiter', ',');
end

function checkIfDir(directory)
    if ~isfolder(directory)
        error('Error. The supplied directory is not a directory');
    end
end