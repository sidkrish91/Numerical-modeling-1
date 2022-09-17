function f = fun_name(t,x)

R = 8.314;
T = 673;
p_ethol = 0.0508; 
p_water = 0.937;
M = p_water / p_ethol;
f = (4.39 * 10^2) * exp ((-2.3 * 10^4)/(R * T)) * (p_ethol ^ 2.42) * ((M - x) ^ 2.71) * ((1 - x) ^ 0.71) * 60;

