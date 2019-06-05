function [rotated_FL, rotated_seg] = rotatePharynx(imStack, segStack, cropDims)
    [centroids, orientations] = calcCentroidsAndOrientations(segStack);
    rotated_FL = imStack;
    rotated_seg = segStack;
    center = size(imStack(:,:,1).') / 2;
    textprogressbar('Translating / Rotating: ');
    for i=1:size(imStack, 3)
        textprogressbar(100 * i / size(imStack, 3));
        translated_FL = imtranslate(rotated_FL(:,:,i), [center(1) - centroids(i, 1) center(2) - centroids(i, 2)]);
        rotated_FL(:,:,i) = imrotate(translated_FL, -orientations(i), 'bilinear', 'crop');
        
        translated_seg = imtranslate(rotated_seg(:,:,i), [center(1) - centroids(i, 1) center(2) - centroids(i, 2)]);
        rotated_seg(:,:,i) = logical(uint8(floor(imrotate(translated_seg, -orientations(i), 'nearest', 'crop'))));
    end
    
    if nargin > 2 
        % cropDims = [height width]
        B = center(2)-cropDims(1)/2;
        T = center(2)+cropDims(1)/2;
        L = center(1)-cropDims(2)/2;
        R = center(1)+cropDims(2)/2;
        
        rotated_FL = rotated_FL(B:T,L:R,:);
        rotated_seg = rotated_seg(B:T,L:R,:);
    end
    
    textprogressbar(' done');
end