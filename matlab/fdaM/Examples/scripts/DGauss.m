function DpG = DGauss(tvec, mu, sigma)

var = sigma.^2;

n = length(tvec);
m = length(mu);
if length(sigma) ~= m
    error('MU and SIGMA not of same length');
end

onesm = ones(1,m);
onesn = ones(n,1);

res   = tvec * onesm - onesn * mu';
expon = res.^2./(2*onesn * var');
DpG   = -res.*exp(-expon);
