function dataSeparated = separateDataByIndexerColumn(data, column, indexer)
% SEPARATEDATABYSTRAIN Given a data matrix (e.g. intensity data) and an
% indexer, this function returns a 1xn cell array (where n is number of
% strains), each cell containing data from only one strain. The order of 
% strains is the same as the order in the indexer.
%
% See also LOADSTANDARDINDEXER
    columnValues = unique(indexer.(column));
    dataSeparated = cell(1, length(columnValues));
    
    for i = 1:length(columnValues)
        columnValue = columnValues(i);
        dataSeparated{i} = data(:, unique(indexer.Animal(indexer.(column) == columnValue)));
    end
end