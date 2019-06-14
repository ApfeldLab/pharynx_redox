function r = genRandInRange(min_,max_,n)
%GENRANDINRANGE Generate n random numbers in range [min_ max].
    r = (max_-min_).*rand(n,1) + min_;
end