function allProps = getMorphProps(segStack)
    allProps = table;
    nAnimals = size(segStack, 3);
    for i=1:nAnimals
        I = segStack(:,:i);
        props = regionprops(I)
    end
end