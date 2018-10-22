function trimmed = trimProfile(intensityProfile)
    trimmed = intensityProfile;
    trimmed(isnan(trimmed))=0;
    for i = 1:size(intensityProfile, 2)
        prof_i = intensityProfile(:,i);

        smoothed = smoothdata(prof_i, 'gaussian', 15);
%         left = find(abs(diff(smoothed)) > 300, 1);
%         right = find(abs(diff(smoothed)) > 200, 1, 'last');
        left = find(smoothed > 2000, 1, 'first');
        right = find(smoothed > 2000, 1, 'last');
        prof_i([1:left right:end]) = 0;
        trimmed(:,i) = prof_i;
    end
    trimmed = ssquare(trimmed);
end