function aggregate = j_aggregate_binary(series)

range=1:size(series,2);
x = 2.^([1:size(series,2)]-1);
aggregate=[];

for n = 1:size(series,1)
    aggregate = cat(1,aggregate,dot(series(n,range),x(range)));
end

    

