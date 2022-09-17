% a script file to solve SK Gupta's example 6.1: Isothermal TRAM (tubular
% reactor with axial mixing) with an irreversible 2nd order reaction 
clc;
clear all;

%bvp syntax
solinit = bvpinit(linspace(0,1,10),[1 0]); %initalization
sol = bvp4c(@rhs_fun_ex6_1, @rhs_bc_ex6_1,solinit); %bvp solver
x = linspace(0,1,7); %intervals to evaluate the differential equations =>mesh size
BS = deval(sol,x); %sol structure evaluator
plot(x,BS(1,:));
xlabel('Axial distance');
ylabel('Concentration');

