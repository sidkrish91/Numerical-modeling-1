clear all;
clc;

beta = 99; %a guess for beta
solinit = bvpinit(linspace(-1,1,10),@bc_init,beta); %solint is is another way to represent init. 
%we can guess zeros as inital guesses.
sol = bvp4c(@rhs_bvp, @bc_bvp, solinit);
x = linspace(-1,1,100); %spliting x into 100 eequal points from 1 to 3
BS = deval(sol,x); %boundary solution to be evaluated from 1 to 3 in 100 equal sized points
plot(x,BS(1,:));

