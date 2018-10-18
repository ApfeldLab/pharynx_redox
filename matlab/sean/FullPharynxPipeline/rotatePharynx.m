function rotated = rotatePharynx(imStack, centroids, orientations)
    rotated = imStack;
    center = size(imStack(:,:,1).') / 2;
    for i=1:size(imStack, 3)
        translated = imtranslate(rotated(:,:,i), [center(1) - centroids(i, 1) center(2) - centroids(i, 2)]);
        rotated(:,:,i) = imrotate(translated, -orientations(i), 'bilinear', 'crop');
    end
end