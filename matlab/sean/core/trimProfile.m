function [trimmed, unscaledBounds, scaledBounds] = trimProfile(intensityProfile, lrBounds)
    trimmed = intensityProfile;
    trimmed(isnan(trimmed))=0;
    
    unscaledBounds = zeros(size(intensityProfile, 2), 2);
    scaledBounds = zeros(size(intensityProfile, 2), 2);
    
    for i = 1:size(intensityProfile, 2)
        prof_i = intensityProfile(:,i);

        smoothed = smoothdata(prof_i, 'gaussian', 15);

        left = find(smoothed >= 2000, 1, 'first');
        right = find(smoothed >= 2000, 1, 'last');
        unscaledBounds(i, :) = [left right];
        
        xs = linspace(lrBounds(i,1), lrBounds(i,2), size(intensityProfile, 1));
        scaledBounds(i, :) = [xs(left) xs(right)];
        
        prof_i([1:left right:end]) = 0;
        trimmed(:,i) = prof_i;
    end
end