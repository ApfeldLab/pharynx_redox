classdef Experiment < handle
    properties
        name
        directory
        metadata
        nAnimals
        strains
        
        raw
        reg
        
        coords410
        coords470
        
        warp
        
        fdConstants
    end
    
    methods
        % Constructor
        function obj = Experiment(directory)
            obj.directory = directory;
            obj.metadata = loadMetadataFromDirectory(directory);
            
            obj.name = obj.metadata.ExperimentID{1};
            
            obj.strains = unique(obj.metadata.Strain);
            obj.nAnimals = size(obj.metadata, 1);
            
            obj.loadRawIntensityData();
            
            obj.coords410 = loadCoords(obj, '410');
            obj.coords470 = loadCoords(obj, '470');
            
            obj.calcRegionMeans('raw', 100);
        end
    end
    
    methods
        function im = getIm410(obj)
            persistent img
            if isempty(img)
                img = loadImage(obj, '410');
            end
            im = img;
        end
        
        function im = getIm470(obj)
            persistent img
            if isempty(img)
                img = loadImage(obj, '470');
            end
            im = img;
        end
        
        function im = getImTL(obj)
            persistent img
            if isempty(img)
                img = loadImage(obj, 'TL');
            end
            im = img;
        end
        
        function registerChannels(obj)
            [obj.reg.fd410, obj.reg.fd470, obj.warp, regInts] =  ...
                ChannelRegister(obj.raw.sq410, obj.raw.sq470, 1000);
            obj.reg.i410 = regInts.m410;
            obj.reg.i470 = regInts.m470;
            
            obj.reg.R = obj.reg.i410 ./ obj.reg.i470;
            obj.reg.OxD = ja_oxd(obj.reg.R);
            obj.reg.E = ja_E(obj.reg.OxD);
            
            obj.calcRegionMeans('reg', 1000);
        end
        
        function idx = filter(obj, boolean_phrase)
            idx = obj.metadata(boolean_phrase,:).Frame;
        end
        
        function plotMeanByStrain(obj, regOrRaw, meas, varargin)
            % TODO: turn off bounded line flag
            % regOrRaw and meas are both strings
            measurement = obj.(regOrRaw).(meas);
            
            groupCol = obj.metadata.Strain;
            groups = unique(groupCol);
            n_groups = numel(groups);
            L = cellfun(@(x)(strcmp(x,groups)), groupCol, 'UniformOutput', false);
            L = cat(2, L{:}).';
            % L is a m x n logical array, where m = # animals, and n = # of groups
            % L(i,j) is 1 if the ith animal is in group j, and 0 otherwise.
            
            resolution = size(measurement, 1);
            
            x = linspace(1,100,resolution);
            means = zeros(resolution, n_groups);
            stds = zeros(resolution, 1, n_groups);
            xs = zeros(resolution, n_groups);
            
            for i=1:numel(groups)
                means(:,i) = mean(measurement(:, L(:,i)).');
                stds(:,1,i) = 1.96*std(measurement(:, L(:,1)),0,2);
                xs(:,i) = x;
            end
            
            if nargin > 3
                ax = varargin{1};
            else
                figure;
                ax = gca();
            end
            
            colors = brewermap(n_groups, 'Set1');
            boundedline(xs, means, stds, 'cmap', colors, 'alpha', 'transparency', 0.15, ax);
            legend(groups);
            
            title(strcat(regOrRaw, meas, ' +/- 1.96*std'));
        end
        
    end
    
    methods (Access = private, Hidden = true)
        function paramStruct = loadParameters(obj)
            param_filepath = fullfile(obj.directory, 'parameters.yaml');
            if ~exist(Name, 'file') == 2
                error('Error. No parameter file in Experiment directory. Please copy parameter file from Image Analysis directory in the shared folder.');
            else 
                paramStruct = ReadYaml(param_filepath);
            end
        end
        
        function loadRawIntensityData(obj)
            obj.raw.i410 = loadIntensity(obj, '410');
            obj.raw.i470 = loadIntensity(obj, '470');
            obj.raw.sq410 = ssquare(clip_sj(obj.raw.i410, 1000));
            obj.raw.sq470 = ssquare(clip_sj(obj.raw.i470, 1000));
            obj.raw.R = obj.raw.sq410 ./ obj.raw.sq470;
            obj.raw.OxD = ja_oxd(obj.raw.R);
            obj.raw.E = ja_E(obj.raw.OxD);
        end
        
        function image = loadImage(obj, suffix)
            fileObj = dir(fullfile(obj.directory, strcat('PA*', suffix, '.tif')));
            image = tiffread2(fullfile(fileObj.folder, fileObj.name));
        end
        
        function md = loadMetadata(obj)
            warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
            mdFile = dir(fullfile(obj.directory, '*.dat'));
            md = readtable(fullfile(obj.directory, mdFile.name), 'Delimiter', ','); % TODO: this is slow (~1sec), optimize at some point
        end
        
        function data = loadIntensity(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '*_intensities.txt')));
            data = dlmread(fullfile(obj.directory, dataFile.name), '', 1, 1);
        end
        
        function coords = loadCoords(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '*_coords.txt')));
            coords = loadCoordinates(fullfile(obj.directory, dataFile.name));
        end
        
        function calcRegionMeans(obj, dataStruct, len)
            % dataStruct` is either 'raw' or 'reg'
            % `len` refers to the sampling size of the data
            %       (the registered data is sampled from the functional data
            %       objects at a variable resolution.)

            obj.(dataStruct).regions.all = table;
            fields = fieldnames(obj.(dataStruct)); % i410, E, ...
            
            for i=1:numel(fields)
                data = obj.(dataStruct).(fields{i});
                if (size(data, 1) == len)
                    % The regions are only defined on the "squared" data
                    % So we only calculate region data on these
                    regions = regionMeans(data);
                    obj.(dataStruct).regions.(fields{i}) = regions;
                    
                    regions.Properties.VariableNames = cellfun(@(x) strcat({x}, '_', fields{i}), regions.Properties.VariableNames);
                    obj.(dataStruct).regions.all = [obj.(dataStruct).regions.all regions];
                end
            end
        end
    end
end