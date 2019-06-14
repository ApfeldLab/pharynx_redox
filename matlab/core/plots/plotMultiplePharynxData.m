function plotMultiplePharynxData(matrices, labels, ylimits, ax, cmap)
%   `matrices` is a cell array of pharynx data
%         pharynx data = PositionInPharynx x Worm
    
    if nargin < 5 || isempty(cmap)
        cmap = cbrewer('qual', 'Set2', max(3, size(matrices, 2)));
    end

    if nargin < 4 || isempty(ax)
        figure;
        ax = gca;
    end
    
%     ts = tinv([0.025  0.975],  size(matrices{1}, 1) - 1);
    
    sz = [size(matrices{1}, 1), sum(cellfun(@(x) size(x, 2), matrices)), size(matrices, 2)]; % X, Y, nMatrices
    
    x = 1:sz(1);
    y_data = zeros(sz(1), sz(3));
    SEM_data = zeros(sz(1), sz(3));
    CI_data = zeros(sz(1), 2, sz(3));
    e_data = zeros(sz(1), 2, sz(3));
    
    interval = 95;
    a = 1 - interval / 100;
    % For each matrix
    for i=1:sz(3)
        L = size(matrices{i}, 1);
        ts = tinv([a/2  1-a/2], L - 1);
        
        y_data(:,i) = mean(matrices{i}, 2);
        SEM_data(:,i) = std(matrices{i}, 0, 2)/sqrt(size(matrices{i}, 2));
        CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
        e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    if nargin < 3 || isempty(ylimits)
        min_ = min(y_data(:), [], 'all');
        max_ = max(y_data(:), [], 'all');
        ylimits = [(min_ - (min_ * .25)) (max_ + (max_ * .25))];
    end
    
    boundedline(x, y_data, e_data, 'cmap', cmap, ax, 'alpha');
    xlim(ax, [1 sz(1)]);
    ylim(ax, ylimits);
    
    if ~isempty(labels)
        legend(labels, 'Interpreter', 'none');
    end
end