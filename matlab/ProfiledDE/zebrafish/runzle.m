function K = runzle(t,XH,Stim)

K = XH;

Tau = 0.8;  %Neural time constants in msec
TauR = 1.9;

K(1) = 1/Tau*(-(17.81 + 47.71*XH(1) + 32.63*XH(1)^2)*(XH(1) - 0.55) - 26*XH(2)*(XH(1) + 0.92) + Stim);  
K(2) = 1/TauR*(-XH(2) + 1.35*XH(1) + 1.03);
