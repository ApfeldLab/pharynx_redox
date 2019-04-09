function midlines = calculateMidlines(TLStack, fluorescenceSegStack, fluorescenceStack, useWeights)
    nAnimals = size(TLStack, 3);
    midlines = cell(nAnimals, 1);
    
    for i=1:nAnimals
%         segTL = imfill(~imbinarize(mat2gray(TLStack(:,:,i))), 'holes');
        segTL = imclose(edge(TLStack(:,:,i),[]), strel('disk',10,6));
        segTL = bwpropfilt(segTL, 'Area', 1);
        eroded = bwmorph(fluorescenceSegStack(:,:,i), 'thin', 2);

        [yTL, xTL] = find(segTL);
        [yFL, xFL] = find(eroded);
        
        stopTLindex = find(xTL > min(xFL), 1, 'first');
        
        x = vertcat(xTL(1:stopTLindex), xFL);
        y = vertcat(yTL(1:stopTLindex), yFL);
        
        im = fluorescenceStack(:,:,i);
        

        ft = fittype('smoothingspline');
        if useWeights
            weights = abs(im(sub2ind(size(im), x, y)));
            opts = fitoptions('Method', 'SmoothingSpline', 'Weights', weights);        
        else
            opts = fitoptions('Method', 'SmoothingSpline');
        end

        opts.SmoothingParam = 0.001;

        [fitresult, ~] = fit(x, y, ft, opts);
        midlines(i) = {fitresult};
    end
end