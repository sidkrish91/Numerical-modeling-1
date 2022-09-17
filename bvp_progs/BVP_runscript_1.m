clear all;
clc;

solinit = bvpinit(linspace(1,3,10),[0 0]); %solint is is another way to represent init. bvpinit picks the inital guess 
%on user input initial guess
sol = bvp4c(@rhs_bvp_1, @bc_bvp_1, solinit);
x = linspace(1,3,100); %spliting x into 100 eequal points from 1 to 3
BS = deval(sol,x); %boundary solution to be evaluated from 1 to 3 in 100 equal sized points
plot(x,BS);

