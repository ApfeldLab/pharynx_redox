function [centroids, orientations] = calcCentroidsAndOrientations(maskStack)
    centroids = zeros(size(maskStack, 3), 2);
    orientations = zeros(size(maskStack, 3), 1);
    maskStack = logical(maskStack);
    for i=1:size(maskStack, 3)
        props = regionprops(maskStack(:,:,i), {'Centroid', 'Orientation'});
        centroids(i, :) = props.Centroid; % 1st coord is x, 2nd is y
        orientations(i) = props.Orientation;
    end
end