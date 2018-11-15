function intensityProfiles = measureIntensityAlongMidlines(fluorescenceStack, midlines, bounds, profileLength, interpMethod)
    nAnimals = size(fluorescenceStack, 3);
    intensityProfiles = zeros(profileLength, nAnimals);
    
    for i=1:nAnimals
        xs = double(bounds(i, 1):bounds(i, 2));
        ys = feval(midlines{i}, xs);
        intensityProfiles(:, i) = improfile(fluorescenceStack(:,:,i), xs, ys, profileLength, interpMethod);
    end
    
    intensityProfiles(isnan(intensityProfiles)) = 0;
end