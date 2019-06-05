function I = loadStandardIndexer(indexerFilePath)
%LOADSTANDARDINDEXER Load the standard (TL/470/410/470/410) indexer
%
%   INPUTS:
%   Takes as input a file path to a STANDARD FORMAT indexer. That is, a csv
%   file with the following structure:
%
%       Strain, Start_Animal, End_Animal
%       HD233,  1,          , 60
%       SAY47,  61,         , 100
%       HD233,  101,        , 200
%
%   The Start/End_Animal columns refer to the ANIMAL and not the absolute
%   image frame (what you would intuitively expect to put in such a file).
%
%   RETURNS:
%   The function returns a table (I for Indexer) with three columsn:
%       ImgFrame, Strain, Lambda, Animal
%
%   ImgFrame gives the absolute index into the interleaved stack
%   Strain refers to the strain of the animal at that frame
%   Lambda refers to the wavelength used at that frame (TL/470_1/410_1/etc)
%   Animal refers to the animal number at that frame
%
%   This structure allows one to use logical operators to index the
%   interleaved stack in an intuitive way.
%   
%   EXAMPLES: 
%   Assume the interleaved stack is called S
%
%   To access all TL images, regardless of strain:
%       S(:, :, I.ImgFrame(I.Lambda == "TL"));
%
%   To access TL images for only HD233:
%       S(:, :, I.ImgFrame(I.Lambda == "TL" & I.Strain == "HD233"));

    indexerTable = readtable(indexerFilePath);
    
    nFramesPerAnimal = 5; % This is true for the "standard" experiment

    nFrames = indexerTable(end,:).End_Animal * nFramesPerAnimal;
    I = table('Size', [nFrames 4], 'VariableTypes', {'uint16', 'string', 'string', 'uint8'}, ...
            'VariableNames', {'ImgFrame', 'Strain', 'Lambda', 'Animal'});

    I.ImgFrame(:) = 1:height(I);
    I.Animal(:) = ceil((1:height(I))/nFramesPerAnimal);
    I.Lambda(1:nFramesPerAnimal:end) = "TL";
    I.Lambda(2:nFramesPerAnimal:end) = "470_1";
    I.Lambda(3:nFramesPerAnimal:end) = "410_1";
    I.Lambda(4:nFramesPerAnimal:end) = "470_2";
    I.Lambda(5:nFramesPerAnimal:end) = "410_2";

    for i=1:height(indexerTable) % for each row in the indexer table
        if i == 1
            start_frame_idx = 1;
        else
            start_frame_idx = (nFramesPerAnimal * indexerTable(i-1,:).End_Animal) + 1;
        end

        end_animal_idx = indexerTable(i,:).End_Animal;
        end_frame_idx = nFramesPerAnimal * end_animal_idx;

        strain = indexerTable.Strain(i,:);

        I.Strain(start_frame_idx:end_frame_idx) = repelem(strain, length(start_frame_idx:end_frame_idx));
    end
end

