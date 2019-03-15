function [distance] = calcMidlineDistance(midlines1, midlines2, scaledBounds1, scaledBounds2, fdObj)
%CALCMIDLINEDISTANCE Summary of this function goes here
%   Detailed explanation goes here

    if length(midlines1) ~= length(midlines2)
        error('Number of Midlines Different');
    end
    prof_len = 100;
    distance = zeros(prof_len, length(midlines1));

    for i=1:length(midlines1)
        xs1 = linspace(scaledBounds1(i, 1), scaledBounds1(i, 2), prof_len);
        xs2 = linspace(scaledBounds2(i, 1), scaledBounds2(i, 2), prof_len);
        ys1 = feval(midlines1{i}, xs1);
        ys2 = feval(midlines2{i}, xs2);

        s = linspace(1,100,prof_len);
        warp_s = eval_fd(s, fdObj(i).warpFD);
        warp_s = warp_s(:,2);

        warp_s_scaled = rescale(warp_s, 1, prof_len);
        warped_ys2 = ys2(min(round(warp_s_scaled), prof_len));
        warped_xs2 = xs2(min(round(warp_s_scaled), prof_len));

        Axy = [xs1; ys1.'];
        Bxy = [warped_xs2; warped_ys2.'];

        distance(:,i) = sqrt(sum((Axy-Bxy).^2,1));
    end
end

