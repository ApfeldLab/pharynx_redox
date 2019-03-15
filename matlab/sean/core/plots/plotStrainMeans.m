function plotStrainMeans(strains, data, title_)
    figure;
    uniqueStrains = unique(strains);
    means = zeros(size(data,1), size(uniqueStrains, 1));
    stds = zeros(size(data,1), 1, size(uniqueStrains, 1));
    
    for i=1:size(uniqueStrains,1)
        strain = uniqueStrains(i);
        dataForStrain = data(:, strains == strain);
        means(:,i) = mean(dataForStrain, 2);
        stds(:, 1, i) = std(dataForStrain, 0, 2);
    end
    
    xs = 1:size(means,1);
    boundedline(xs, means, stds, 'alpha');
    legend(cellstr(uniqueStrains));
    title(title_);
end