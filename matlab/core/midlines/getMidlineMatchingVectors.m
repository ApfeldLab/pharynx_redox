function [vecs] = getMidlineMatchingVectors(midlines1, midlines2, scaledBounds1, scaledBounds2, fdObj, nVecs)
    if length(midlines1) ~= length(midlines2)
        error('Number of Midlines Different');
    end
   
    vecs = cell(length(midlines1), 1);

    for i=1:length(midlines1)
        xs1 = linspace(scaledBounds1(i, 1), scaledBounds1(i, 2), nVecs);
        xs2 = linspace(scaledBounds2(i, 1), scaledBounds2(i, 2), nVecs);
        ys1 = feval(midlines1{i}, xs1);
        ys2 = feval(midlines2{i}, xs2);

        s = linspace(1,100,nVecs);
        warp_s = eval_fd(s, fdObj(i).warpFD);
        warp_s = warp_s(:,2);

        warp_s_scaled = rescale(warp_s, 1, nVecs);
        warped_ys2 = ys2(min(round(warp_s_scaled), nVecs));
        warped_xs2 = xs2(min(round(warp_s_scaled), nVecs));
        
%         theta = atan((warped_ys2 - ys1) ./ (warped_xs2.' - xs1.'));
%         u = cos(theta);
%         v = sin(theta);
        vecs{i} = [(warped_xs2.' - xs1.') (warped_ys2 - ys1)];
    end
end
