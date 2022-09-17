% to solve 1 Dimensional Steady state Pseudo homogeneous model

function dz=SS_pseudohomo_batch_1(z,x)

%declaring constants
 Us = 1; %superficial velocity
 L = 3; % tube length
 Dt = 2.54*1.e-2; %tube diameter
 dp = 3 * 1.e-3; %particle diameter
 row_b = 1300;
 M = 29.48;
 delh = 1285409 * 1.e+3;
 Cp = 0.992 * 1.e+3;
 Ptot = 1;
 Pbo = 0.211;
 U = 0.096 * 1.e+;
 Ta = 352 + 273;
 row_g = Ptot * M / (0.08206 * Ta);
 