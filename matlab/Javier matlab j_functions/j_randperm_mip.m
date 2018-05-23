function a = j_randperm_mip(num, series1, series2)
%function a = j_randperm_mip(num, series1, series2)
%generates a column vector of size num where each element is
%the j_mutualinfopair of series1 and series2, 
%(with series2 in random order)

a = [];
for n = 1:num
    index = randperm(size(series2,1));
    a(n,1) = j_mutualinfopair(series1, series2(index));
end

    