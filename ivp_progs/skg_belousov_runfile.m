%to solve the belousov ode
%initial values
clear;clc;
yo = [4,1.1,4];
tspan = [0 5];
options=odeset('JPattern','on');
[t,y]=ode15s(@skg_belousov_rxns,tspan,yo,options);
plot(t,y);