function plotPharynxDataByMovement(r, mvmt_mode, ax)
    x=1:100;
    ts = tinv([0.025  0.975],size(r, 1)-1);

    y_data = zeros(100, 3);
    SEM_data = zeros(100, 3);
    CI_data = zeros(100, 2, 3);
    e_data = zeros(100, 2, 3);
    
    for i=1:4
        idx = mvmt_mode == i - 1;
        y_data(:,i) = mean(r(:,idx), 2).';
        SEM_data(:,i) = std(r(:,idx), 0, 2)/sqrt(size(r, 1));
        CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
        e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    
    cmap_scaled = flipud(cbrewer('qual', 'Dark2', 4));
    if nargin > 2
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled, ax);
        xlim(ax, [1 100]);
        ylim auto;
    else
        figure;
        boundedline(x, y_data, e_data, 'cmap', cmap_scaled);
        xlim([1 100]);
        ylim auto;
    end
    
    legend('0', '1', '2', '3', 'Location', 'northwest');
    xlabel('Pharynx Space');
end