function plotMultiplePharynxData(matrices, labels, ylimits, ax)
%   `matrices` is a cell array of pharynx data
%         pharynx data = PositionInPharynx x Worm
    
%     matrices = horzcat(matrices{:});

    if nargin < 4 || isempty(ax)
        figure;
        ax = gca;
    end
    
%     if nargin < 3 || isempty(ylimits)
%         ylimits = [(min(matrices(:)) - (min(matrices(:) * .25))) (max(matrices(:)) + (max(matrices(:) * .25)))];
%     end
    
    ts = tinv([0.025  0.975],  size(matrices{1}, 1) - 1);
    
    sz = [size(matrices{1}, 1), sum(cellfun(@(x) size(x, 2), matrices)), size(matrices, 2)]; % X, Y, nMatrices
    
    x = 1:sz(1);
    y_data = zeros(sz(1), sz(3));
    SEM_data = zeros(sz(1), sz(3));
    CI_data = zeros(sz(1), 2, sz(3));
    e_data = zeros(sz(1), 2, sz(3));
    
    % For each matrix
    for i=1:sz(3)
            y_data(:,i) = mean(matrices{i}, 2);
            SEM_data(:,i) = std(matrices{i}, 0, 2)/sqrt(size(matrices{i}, 2));
            CI_data(:,:,i) = y_data(:,i) + ts.*SEM_data(:,i);
            e_data(:,:,i) = abs(y_data(:,i) - CI_data(:,:,i));
    end
    
    cmap = cbrewer('qual', 'Dark2', max(3, size(matrices, 3))); 
    
    boundedline(x, y_data, e_data, 'cmap', cmap, ax, 'alpha');
    xlim(ax, [1 sz(1)]);
    ylim(ax, ylimits);
    
    if ~isempty(labels)
        legend(labels);
    end
end