function seg = segmentPharynx(imStack, useAC, thresh)
    seg = zeros(size(imStack));

    acSE = strel('square', 3);

    steps = size(imStack, 3);
    disp('Segmenting Images... Please hold');
    h = fspecial('unsharp');
    for i=1:steps
        mask = imfill(imdilate(imclose(edge(imStack(:,:,i), 'sobel'), strel('diamond',2)), strel('diamond',2)), 'holes');
        mask = bwpropfilt(logical(mask), 'Area', 1);
        
        if useAC
            unsharped = imfilter(imStack(:,:,i), h);
            seg(:,:,i) = imdilate(activecontour(unsharped, mask, 15, 'Chan-Vese', 'SmoothFactor', 3), acSE);
            seg(:,:,i) = bwpropfilt(logical(seg(:,:,i)), 'Area', 1);
        else
            seg(:,:,i) = mask;
        end
        
        if nargin > 2
            mask(mask < thresh) = 0;
        end
    end
    seg = logical(seg);
    disp('DONE Segmenting Images');
end