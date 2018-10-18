function intensityProfiles = measureIntensityAlongMidlines(fluorescenceStack, midlines)
    nAnimals = size(fluorescenceStack, 3);
    imWidth = size(fluorescenceStack, 2);
    intensityProfiles = zeros(imWidth, nAnimals);
    xs = 1:imWidth;
    
    for i=1:nAnimals
        intensityProfiles(:, i) = improfile(fluorescenceStack(:,:,i), xs, feval(midlines{i}, xs), imWidth);
    end
%     intensityProfiles = ssquare(clip_sj(intensityProfiles, threshold));
end