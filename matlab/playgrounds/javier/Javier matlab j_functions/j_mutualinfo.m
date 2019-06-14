function I = j_mutualinfo(series)

I = zeros(size(series,2));

for n = 1:size(series,2);
    for m = 1:size(series,2);
        I(n,m) = j_mutualinfopair(series(:,n),series(:,m));
    end
end

%I

 