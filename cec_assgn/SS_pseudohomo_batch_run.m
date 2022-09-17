% program to run SS_pseudohomo_batch

clc;
clear;

z = [0 3];
x0 = [0.015 625];
[z x] = ode15s(@SS_pseudohomo_batch,z,x0);

 