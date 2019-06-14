function trimmed = trimProfile(intensityProfile, threshold)
    trimmed = intensityProfile;
    trimmed(isnan(trimmed))=0;
    
    for i = 1:size(intensityProfile, 2)
        prof_i = intensityProfile(:,i);

        smoothed = smoothdata(prof_i, 'gaussian', 15);

        left = find(smoothed >= threshold, 1, 'first');
        right = find(smoothed >= threshold, 1, 'last');
        
        prof_i([1:left right:end]) = 0;
        trimmed(:,i) = prof_i;
    end
end