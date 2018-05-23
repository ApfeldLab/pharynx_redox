function [a, names] = j_make_strain_cellarray(series1, series2)
%function [a, names] = j_make_strain_cellarray(series1, series2)
%makes a cell array of tables (row = worm#, col1=time to event, col2=
%censoring) one for each element of series2.
%names is a column vector with the names of the elements of series2 from
%which each cell is derived

sortrows(series2,1);
a{1}= [];
names = [];
for n = 1:size(series2,1)
    index = find(names == series2(n,1));
    if isempty(index) == 1
        names = [names; series2(n,1)];
        index2 = find(series2(:,1) == series2(n,1));
        a{length(names)}= series1(index2,:);
    end
end
