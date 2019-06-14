function a = j_series_randperm(series)
%function a = j_series_randperm(series)
%performs random permutation of each column

a = [];
for n = 1:size(series,2)
    index = randperm(size(series,1));
    a(:,n) = series(index,n);
end
