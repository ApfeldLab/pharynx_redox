classdef Experiment < handle
    methods
        % Constructor
        function obj = Experiment(dirPath)
            % TODO:
            % Check if a .mat file exists in the folder. If it does, just
            % load that .mat file as the class object.
            
            if (~endsWith(dirPath, filesep))
                dirPath = strcat(dirPath, filesep);
            end
            obj.varname = evalin('caller','inputname(1)');
            
            obj.directory = dirPath;
            
            dirparts = split(dirPath, filesep);
            obj.name = dirparts{end-1};
            obj.metadata = loadMetadata(obj);
            obj.strains = unique(obj.metadata.Strain);
            obj.nAnimals = size(obj.metadata, 1);
            
            obj.raw.i410 = loadIntensity(obj, '410');
            obj.raw.i470 = loadIntensity(obj, '470');
            obj.raw.sq410 = ssquare(clip_sj(obj.raw.i410, 1000));
            obj.raw.sq470 = ssquare(clip_sj(obj.raw.i470, 1000));
            obj.raw.R = obj.raw.sq410 ./ obj.raw.sq470;
            obj.raw.OxD = ja_oxd(obj.raw.R);
            obj.raw.E = ja_E(obj.raw.OxD);
            
            obj.coords410 = loadCoords(obj, '410');
            obj.coords470 = loadCoords(obj, '470');
            
            obj.calcRegionMeans('raw', 100);
            
            obj.fdConstants = struct(...
                'smoothL', 0.0891, 'smoothOrder', 6, 'smoothBasis', 96, ...
                'warpL',     5000, 'warpOrder',   4,   'warpBasis',   6 ...
                );
        end
    end
    
    properties
        varname
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
        function im = getIm410(obj)
            persistent img
            if isempty(img)
                img = loadImage(obj, '410_subMode');
            end
            im = img;
        end
        
        function im = getIm470(obj)
            persistent img
            if isempty(img)
                img = loadImage(obj, '470_subMode');
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
                ChannelRegister(obj.raw.sq410, obj.raw.sq470);
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
    end
    
    methods (Access = private, Hidden = true)
        function image = loadImage(obj, suffix)
            fileObj = dir(fullfile(obj.directory, strcat('PA*', suffix, '.tif')));
            image = tiffread2(fullfile(fileObj.folder, fileObj.name));
        end
        
        function md = loadMetadata(obj)
            warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames'); % TODO do i need to worry about this warning? Suppresed for now
            mdFile = dir(fullfile(obj.directory, '*.dat'));
            md = readtable(fullfile(obj.directory, mdFile.name)); % TODO: this is slow (~1sec), optimize at some point
        end
        
        function data = loadIntensity(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '_subMed_intensities.txt')));
            data = dlmread(fullfile(obj.directory, dataFile.name), '', 1, 1);
        end
        
        function coords = loadCoords(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '_subMed_coords.txt')));
            coords = loadCoordinates(fullfile(obj.directory, dataFile.name));
        end
        
        function calcRegionMeans(obj, dataStruct, len)
            % 1dataStruct` is either 'raw' or 'reg'
            % `len` refers to the sampling size of the data
            %       the registered data is sampled from the functional data
            %       objects at a variable resolution.

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