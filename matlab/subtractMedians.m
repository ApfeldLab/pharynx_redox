function subMedians = subtractMedians(imgStack)
    medians = median(imgStack, [1 2]);
    subMedians = zeros(size(imgStack), 'double');
    for i=1:size(imgStack, 3)
        img = imgStack(:,:,i);
        subMedians(:,:,i) = img - medians(i);
    end
end