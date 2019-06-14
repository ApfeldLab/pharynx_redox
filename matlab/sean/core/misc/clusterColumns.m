% When I image on MetaMorph, sometimes I accidentally zoom past a worm. If
% I used the # worms per column that I recorded when plating the animals,
% then the mapping from column to strain might be messed up.

% So what this function does is take in the coordinates that you can get 
% from running the log-coordinates journal on MetaMorph, cluster them into
% however many columns you specify, then spit out the column # for each
% animal.

% It also plots the coordinates where color is mapped to column #
% as a sanity check. If you see a column with more than 1 color, then
% something went wrong in the cluster and you should ignore the output of
% this function and just try to assign the columns manually.

% `reindexed_cols` is populated where each row is given the value of the
% column it belongs to. This can be saved as a CSV file as such:
%   csvwrite('path/to/filename.csv', reindexed_cols)
% This CSV can then be pasted into the "Column" column on JMP.

% TODO: return column #s as well
function strains = clusterColumns(camera_position_filename, columnStrains, nChannels, shouldPlot)
    t = readtable(camera_position_filename);
    X = t.PlaneStagePositionX(1:nChannels:end);
    
    cols = kmeans(X, numel(columnStrains));
    clusters = unique(cols, 'stable'); % stable means don't sort

    strains = cell(numel(cols),1);
    for i=1:max(clusters)
        strains(cols==clusters(i)) = {columnStrains{i}};
    end

    if shouldPlot
        [~, ~, group_idx] = unique(strains);
        Y = t.PlaneStagePositionY(1:nChannels:end);
        scatter(X, Y, 30, group_idx);
    end
end