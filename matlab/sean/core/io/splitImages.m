function [im410_1, im410_2, imTL] = splitImages(allImages, group_id_str, indexer)
    % NOTE: The way these images are split are specific to this experiment!
    % Read the documentation on the loadIndexer method for more info on
    % splitting up images with the indexer.
    %
    % For this particular analysis, we are interested in the 410/410 data set.
    % I am also splitting up this data on a per-animal basis.
    im410_1 = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "410/410" & indexer.SetFrame == 1 & indexer.Strain == group_id_str));
    im410_2 = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "410/410" & indexer.SetFrame == 2 & indexer.Strain == group_id_str));

    imTL = allImages(:,:,indexer.ImgFrame(indexer.LambdaGroup == "TL" & indexer.Strain == group_id_str));
    imTL = repmat(imTL, [1 1 size(im410_1, 3)]);
end