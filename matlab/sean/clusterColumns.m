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

function reindexed_cols = clusterColumns(filename, nCols, shouldPlot)
    
    t = readtable(filename);
    X = t.PlaneStagePositionX;
    cols = kmeans(X, nCols);
    clusters = unique(cols, 'stable'); % don't want unique to sort this

    reindexed_cols = cols;
    for i=1:max(clusters)
        reindexed_cols(cols==clusters(i)) = i;
    end

    if shouldPlot
        scatter(t.PlaneStagePositionX, t.PlaneStagePositionY, 30, reindexed_cols);
    end
    
end