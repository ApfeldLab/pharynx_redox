function intensityProfiles = measureIntensityAlongMidlines(...
    fluorescenceStack, midlines, profileLength, interpMethod)

    imWidth = size(fluorescenceStack, 2);
    nAnimals = size(fluorescenceStack, 3);
    intensityProfiles = zeros(profileLength, nAnimals);
    
    bounds_ = [1 imWidth];
    
    textprogressbar('measuring under midlines: ');
    for i=1:nAnimals
        textprogressbar(100 * (i/nAnimals));
        xs = double(bounds_(1):bounds_(2)).';
        ys = feval(midlines{i}, xs);
        intensityProfiles(:, i) = improfile(fluorescenceStack(:,:,i), xs, ys, profileLength, interpMethod);
    end
    textprogressbar('done');
    
    intensityProfiles(isnan(intensityProfiles)) = 0;
end