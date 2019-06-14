lnorm = square(data);
kym = getint(lnorm);

lnorm = square(data);
kym = getint(lnorm);
R= gdivide (kym_410, kym_470);
R_mean= mean(R, 1);
R_plusmean = vertcat(R_mean, R)

% R_mean= mean(R_WT_10do, 1);
% R_plusmean = vertcat(R_mean, R_WT_10do);