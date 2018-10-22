function intensityProfiles = measureIntensityAlongMidlines(fluorescenceStack, midlines)
    nAnimals = size(fluorescenceStack, 3);
    imWidth = size(fluorescenceStack, 2);
    intensityProfiles = zeros(imWidth * 10, nAnimals);
    xs = 1:imWidth;
    
    for i=1:nAnimals
        ys = feval(midlines{i}, xs);
        intensityProfiles(:, i) = improfile(fluorescenceStack(:,:,i), xs, ys, imWidth * 10, 'bilinear');
    end
    intensityProfiles(isnan(intensityProfiles)) = 0;
%     intensityProfiles = ssquare(clip_sj(intensityProfiles, threshold));
end