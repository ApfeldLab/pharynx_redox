function dataSeparated = separateDataByStrain(data, indexer)
% SEPARATEDATABYSTRAIN Given a data matrix (e.g. intensity data) and an
% indexer, this function returns a 1xn cell array (where n is number of
% strains), each cell containing data from only one strain. The order of 
% strains is the same as the order in the indexer.
%
% See also LOADSTANDARDINDEXER
    strains = unique(indexer.Strain);
    dataSeparated = cell(1, length(strains));
    
    for i = 1:length(strains)
        strain = strains(i);
        dataSeparated{i} = data(:, unique(indexer.Animal(indexer.Strain == strain)));
    end
end