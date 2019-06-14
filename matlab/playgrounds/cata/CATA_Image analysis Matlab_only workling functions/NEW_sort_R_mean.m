% R_daf2_2do= gdivide (kym_daf2_2do_410, kym_daf2_2do_470);
% 
% R_mean= mean(R_daf2_2do, 1);
% 
% R_plusmean = vertcat(R_mean, R_daf2_2do);

%R= gdivide (kym_TB_410, kym_TB_470);

R= gdivide (kym_410, kym_470);

R_mean= mean(R, 1);

R_plusmean = vertcat(R_mean, R);