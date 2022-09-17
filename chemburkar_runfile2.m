%% Change in gamma(5) value to check for parametric
%%%%%%%%%%%%%%sensitivity%%%%%%%%
clc;clear;
close all;
% global gamma;

y0 = [0.0213 0.0375 4.629]; %initial values for dimensionless variables
tspan = 0:1e-4:25; %step size of 10^-4. Ref: Chemburkar et al for step size lowering for period doubling

gamma = [0.26 1e307 0.50 1.0 19.85565 0.00 57.77 -0.426];
[t2,y2] = ode23t(@(t,y) chemburkar_func(t,y,gamma),tspan,y0);
%% Plotting of curves
f2=figure();
f2=plot(t2,y2(1:end,3));
title('Temperature(dimensionless) vs Time');
legend (f2,'gamma(5) = 19.85565');  