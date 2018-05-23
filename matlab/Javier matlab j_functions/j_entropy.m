function H = j_entropy(series)

%output is the entropy of each column of the input series

%series=[1 2 3 4; 1 6 7 8]';
%series = [6 5 4 3; 1 2 3 5]'


%find all entropies
H = zeros(1,size(series,2));
for n = 1:size(series,2)
    cat_count = j_count(series(:,n));
    p = cat_count(:,2)./sum(cat_count(:,2));
    H(n) = -sum(log2(p.^p));
end

%H



    
    