classdef (Abstract) Experiment
    properties
        name
        experimenter
        imagingDate
        indexer
        
        rawImageStack
    end
    
    properties (Constant)
        requiredVersion = '9.5.0.944444 (R2018b)';
    end
    
    methods
       % Constructor
       function E = Experiment(imageFilePath, indexerTable, experimenter, imagingDate, name)
           checkVersion();
           
           E.name = name;
           E.imagingDate = imagingDate;
           E.experimenter = experimenter;
           E.indexer = loadIndexer(indexerTable);
           E.rawImageStack = E.loadImageStack(imageFilePath);
       end
    end
    
    methods (Abstract)
        loadIndexer(indexerTable)
    end
    
    methods (Static)
        function checkVersion
            % Ensure that MATLAB is the correct version. Currently, we
            % force the version to be the EXACT same as the development
            % version. This is because MathWorks sometimes introduces API
            % changes that break the code with versions.
            assert(strcmp(version, Experiment.requiredVersion), ...
                sprintf('Incorrect MATLAB version. Must be [%s]. Currently running verion [%s].', Experiment.requiredVersion, version))
        end
    end
    
    methods (Access = private)
        function S = loadImageStack(tiffImageStackFilePath)
            % Read a multi-frame TIFF image stack into a uint16 matrix
            % (width x height x nFrames) 
    
            InfoImage=imfinfo(tiffImageStackFilePath);
            mImage=InfoImage(1).Width;
            nImage=InfoImage(1).Height;
            NumberImages=length(InfoImage);
            S=zeros(nImage, mImage, NumberImages, 'uint16');

            TifLink = Tiff(tiffImageStackFilePath, 'r');
            for i=1:NumberImages
               TifLink.setDirectory(i);
               S(:,:,i)=TifLink.read();
            end
            TifLink.close();
        end
    end
end