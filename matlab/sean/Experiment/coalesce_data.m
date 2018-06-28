DIRECTORY = '/Users/sean/Desktop/2018_06_14_SAY98_HD233';

experimentID = getExperimentID(DIRECTORY);
measurementTable = loadMeasurements(DIRECTORY);
movementRep = loadMvmtRep(DIRECTORY);

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

function mvmtRep = loadMvmtRep(directory)
    mvmtRep = loadDataFile(directory, '*_mvmt.csv');
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