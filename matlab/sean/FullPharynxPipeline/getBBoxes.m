function bBoxes = getBBoxes(segStack)
    nAnimals = size(segStack, 3);
    bBoxes = zeros(nAnimals, 4);
    for i=1:nAnimals
        I = segStack(:,:,i);
        stats = regionprops(I, 'BoundingBox');
        bBoxes(i,:) = stats.BoundingBox;
    end
end