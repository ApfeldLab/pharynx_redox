function I = loadIndexer(indexerFilePath)
%LOADINDEXER Load the indexer from the given filepath
%
%   The indexer allows you to index frames of a TIFF stack according to
%   metadata stored in a .csv file. The .csv file should be structured as
%   such:
%       start:end,group,lambda1/lambda2/lambda3 (etc)
%   
%   So, for example, the file might look like this:
%       1:1,     Animal 1, TL
%       2:127,   Animal 1, 470/410
%       128:241, Animal 1, 410/410
%       242:242, Animal 2, TL
%       243:364, Animal 2, 470/410
%       365:488, Animal 2, 410/410
%
%   Or, maybe like this:
%       1:90,   HD233, TL/470/410
%       91:211, SAY47, TL/410/410
%
%   The start:end column MUST contain a start AND end, even if it is a
%   single frame. In that case, use the same frame for start and end (see
%   example 1, rows 1 and 4).
%
%   -----------------------------------------------------------------------
%   RETURNS:
%   A table with as many rows as the largest `end` in the indexer file, and
%   the following columns:
%       ImgFrame, Strain, LambdaGroup, Lambda, SetFrame
%   
%   ImgFrame [uint8]:
%       The index of the TIFF stack
%
%   Strain [string]:
%       The group from the indexer file (Animal 1, Animal 2, HD233, etc.)
%
%   LambdaGroup [string]:
%       The group of lambdas that you took your images at. This simply
%       mirrors the lambda1/lambda2/... that you put in the indexer file
%
%   Lambda [string]:
%       The actual lambda of the row
%
%   SetFrame [uint8]:
%       Which image in the burst this row corresponds to (so if you took a
%       3 channel burst such as TL/410/410, this column lets you
%       distinguish the first and second 410nm frames)
%
%   -----------------------------------------------------------------------
%   USAGE:
%   Say you have a TIFF stack that just came from metamorph. This indexer
%   allows you to select particular frames from that stack according to
%   filters that you can build.
%
%   For example, suppose I just imaged a bunch of animals from two
%   genotypes using TL/470/410/470/410. Let the name of the object that this
%   function returns be `I`. Let the TIFF stack be `stack`.
%   
%   Access all TL images:
%   all_TL = stack(:, :, I.ImgFrame(I.Lambda == "TL"));
%
%   Access all TL images from a particular genotype:
%   all_TL_HD233 = stack(:, :, I.ImgFrame(I.Lambda == "TL" & I.Strain == "HD233"));
%   all_TL_SAY47 = stack(:, :, I.ImgFrame(I.Lambda == "TL" & I.Strain == "SAY47"));
%
%   Access all of the first 470 and first 410 images from all genotypes:
%   im470 = stack(:,:, I.ImgFrame(I.SetFrame == 2);
%   im410 = stack(:,:, I.ImgFrame(I.SetFrame == 3);
%
%   -----------------------------------------------------------------------
%   Try playing around with different queries to get a feel for how to
%   effectively access the data!
       
    indexer = readtable(indexerFilePath, 'Delimiter', ',', 'Format', '%s%s%s', 'ReadVariableNames', false);

    % This line is ugly, but essentially looks down the first column of the
    % index file, converts all the strings to numbers, and finds the
    % maximum
    nFrames = max(cell2mat(cellfun(@(x) cellfun(@str2double, strsplit(x, ':')), indexer.Var1, 'UniformOutput', false)), [], 'all');
    
    I = table('Size', [nFrames 5], 'VariableTypes', {'uint16', 'string', 'string', 'string', 'uint8'}, ...
        'VariableNames', {'ImgFrame', 'Strain', 'LambdaGroup', 'Lambda', 'SetFrame'});

    for i=1:height(indexer)
        channels = strsplit(indexer.Var3{i}, '/');
        bounds = uint16(cellfun(@str2double, strsplit(indexer.Var1{i}, ':')));

        I.Strain(bounds(1):bounds(2)) = indexer.Var2{i};
        I.LambdaGroup(bounds(1):bounds(2)) = indexer.Var3{i};
        I.ImgFrame(bounds(1):bounds(2)) = bounds(1):bounds(2);

        nChannels = length(channels);
        for j=1:nChannels
            I.Lambda(bounds(1)+j-1:nChannels:bounds(2)) = channels(j);
            I.SetFrame(bounds(1)+j-1:nChannels:bounds(2)) = j;
        end
    end
end

