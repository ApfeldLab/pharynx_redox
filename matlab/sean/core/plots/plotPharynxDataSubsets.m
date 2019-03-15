function plotPharynxDataSubsets(data3d, ax)
%     Each group should be a new 2d plane of pharynx data, 
%       (Worms x Pharynx x Group)
    x=1:size(data3d, 1);
    ts = tinv([0.025  0.975],  size(data3d, 1) - 1);

    y_data = zeros(100, 3);
    SEM_data = zeros(100, 3);
    CI_data = zeros(100, 2, 3);
    e_data = zeros(100, 2, 3);
    
    for i=1:size(data3d, 3)
        y_data(:,i) = mean(data3d(:,:,i), 2).';
        SEM_data(:,i) = std(data3d(:,:,i), 0, 2)/sqrt(size(data3d, 1));
        CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
        e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    cmap_ = flipud(cbrewer('qual', 'Dark2', max(3, size(data3d, 3))));
    
    if nargin > 1
        boundedline(x, y_data, e_data, 'cmap', cmap_, ax);
        xlim(ax, [1 100]);
    else
        figure;
        boundedline(x, y_data, e_data, 'cmap', cmap_);
        xlim([1 100]);
    end 
end