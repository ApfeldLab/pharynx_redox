function midlines = calculateMidlines(TLStack, fluorescenceMaskStack)
    nAnimals = size(TLStack, 3);
    midlines = cell(nAnimals, 1);
    
    for i=1:nAnimals
        segTL = imfill(~imbinarize(mat2gray(TLStack(:,:,i))), 'holes');
        segTL = bwpropfilt(segTL, 'Area', 1);
        eroded = bwmorph(fluorescenceMaskStack(:,:,i), 'thin', 2);

        [yTL, xTL] = find(segTL);
        [yFL, xFL] = find(eroded);
        
        stopTLindex = find(xTL > min(xFL), 1, 'first');
        
        x = vertcat(xTL(1:stopTLindex), xFL);
        y = vertcat(yTL(1:stopTLindex), yFL);

        ft = fittype('smoothingspline');
        opts = fitoptions('Method', 'SmoothingSpline');
        opts.SmoothingParam = 0.001;

        [fitresult, ~] = fit(x, y, ft, opts);
        midlines(i) = {fitresult};
    end
end