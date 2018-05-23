function D = j_infodist1(series)

% D = 1 - I(a,b)./(Ha U Hb) | a,b are columns of the input series

%series=[1 2 3 4; 1 6 7 8]';
%series = [6 5 4 3; 1 2 3 5]'


%find all entropies and mutual information pairs
H = j_entropy(series);
I = j_mutualinfo(series);

%calculate distance
D = zeros(size(I));
for n = 1:size(D,1)
    for m = 1:size(D,2)
        D(n,m)=I(n,m)/(H(n)+H(m)-I(n,m));
    end
end
D= 1-D;
        
