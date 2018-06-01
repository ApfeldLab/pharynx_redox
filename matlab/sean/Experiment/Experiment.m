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
        end
    end
    
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
            [obj.reg.fd410, obj.reg.fd470, obj.warp, regInts] = ChannelRegister(obj.raw.sq410, obj.raw.sq470);
            obj.reg.i410 = regInts.m410;
            obj.reg.i470 = regInts.m470;
            
            obj.reg.R = obj.reg.i410 ./ obj.reg.i470;
            obj.reg.OxD = ja_oxd(obj.reg.R);
            obj.reg.E = ja_E(obj.reg.OxD);
        end
    end
    
    methods (Access = private, Hidden = true)
        function image = loadImage(obj, suffix)
            fileObj = dir(fullfile(obj.directory, strcat('PA*', suffix, '.tif')));
            image = tiffread2(fullfile(fileObj.folder, fileObj.name));
        end
        
        function md = loadMetadata(obj)
            mdFile = dir(fullfile(obj.directory, '*.dat'));
            md = readtable(fullfile(obj.directory, mdFile.name));
        end
        
        function data = loadIntensity(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '_intensities.txt')));
            data = dlmread(fullfile(obj.directory, dataFile.name));
        end
        
        function coords = loadCoords(obj, channel)
            dataFile = dir(fullfile(obj.directory, strcat('*', channel, '_coords.txt')));
            coords = loadCoordinates(fullfile(obj.directory, dataFile.name));
        end
    end
end