function plotMultipleMidlines(ratioImage, flImage, mask, midline1, midline2, ...
    scaledBounds1, scaledBounds2, warpFD)
    
%     I = ratioImage;
%     I_masked = I.*segImage;
%     I_masked(I_masked == 0) = NaN;
%     imAlpha=ones(size(I_masked));
%     imAlpha(isnan(I_masked))=0;
% 
% %     imagesc(ax, I, 'AlphaData', imAlpha); 
%     imagesc(ax, I);
%     colormap(cmap); colorbar;
%     caxis([.75 1.25]);
%     set(gca,'color',0*[1 1 1]);


    cmap_ = cbrewer('div', 'RdBu', 256, 'PCHIP');
    R_adj = ja_adjust_brightness(ratioImage, flImage, 1200, cmap_, .85, 1.25);
    imshow(R_adj);
    
    axis square;
    
    bbox = regionprops(mask, 'BoundingBox');
    bbox = floor(bbox.BoundingBox);
    buffer = 15;

    left = bbox(1) - buffer;
    top = bbox(2) - buffer;
    right = bbox(1) + bbox(3) + buffer;
    bottom = bbox(2) + bbox(4) + buffer;
    xlim([left right]);
    ylim([top bottom]);
    
    hold on;
    plot(midline1, 'r-');
    plot(midline2, 'g-');
    
    % Concommitant Points
    prof_len = 100;
    xs1 = linspace(scaledBounds1(1), scaledBounds1(2), prof_len);
    xs2 = linspace(scaledBounds2(1), scaledBounds2(2), prof_len);
    
    ys1 = feval(midline1, xs1);
    ys2 = feval(midline2, xs2);
    
    s = linspace(1,100,prof_len);
    warp_s = eval_fd(s, warpFD);
    warp_s = warp_s(:,2);

    warp_s_scaled = rescale(warp_s, 1, prof_len);
    warped_ys2 = ys2(min(round(warp_s_scaled), prof_len));
    warped_xs2 = xs2(min(round(warp_s_scaled), prof_len));
    
    for i=1:length(xs1)
        plot([xs1(i) warped_xs2(i)], [ys1(i) warped_ys2(i)], 'k-');
    end
    
    lgnd = legend('410nm Midline', '470nm Midline');
    set(lgnd, 'color', 'w');
    
    hold off;
end