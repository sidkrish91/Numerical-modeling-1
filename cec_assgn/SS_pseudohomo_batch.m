% to solve 1 Dimensional Steady state Pseudo homogeneous model

function dy=SS_pseudohomo_batch(z,x)

dy = zeros(2,1);

Pa = x(1);
T = x(2);

%declaring constants

 Us = 1; %superficial velocity
 L = 3; % tube length
 Dt = 2.54*1.e-2; %tube diameter
 dp = 3 * 1.e-3; %particle diameter
 row_b = 1300;
 M = 29.48;
 delH = 1285409 * 1.e+3;
 Cp = 0.992 * 1.e+3;
 Ptot = 1;
 Pbo = 0.211;
 U = 0.096 * 1.e+3;
 Ta = 625;
 row_g = Ptot * M / (0.08206 * Ta);
 k = exp(19.837 - 13696/T);
 
 
 dy(1) = k * row_b * Pa * Pbo * 0.08206 * T;
 dy(2) = (delH * row_b * k * Pa * Pbo -  (4 * U *(T - Ta) / Dt)) / (Us * row_g);
 
 