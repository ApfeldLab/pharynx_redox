function plotMovmentStratification(movementLabels, data, title_)
    figure;
    
    uniqueLabels = unique(movementLabels);
    means = zeros(size(data,1), size(uniqueLabels, 1));
    stds = zeros(size(data,1), 1, size(uniqueLabels, 1));
    
    for i=1:size(uniqueLabels,1)
        label = uniqueLabels(i);
        dataForLabel = data(:, movementLabels == label);
        means(:,i) = mean(dataForLabel, 2);
        stds(:, 1, i) = std(dataForLabel, 0, 2);
    end
    
    xs = 1:size(means,1);
    boundedline(xs, means, stds, 'alpha');
    legend(cellstr(num2str(uniqueLabels)));
    title(title_);
end