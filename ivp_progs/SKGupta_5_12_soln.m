%solution file for SKGupta_5_12_plug
clc
clear all;
close all;

%intializing input variables
y0 = [0.002;19.1071;0;0;0;0;0;0;50];
tspan = 0:0.001:35;

%running ODE15s
[t,y] = ode15s(@SKGupta_5_12_plug, tspan, y0,10e-7);

%plotting curves

plot(12.81 * t,y(:,9));
hold on
[yax,y9,mu]=plotyy(12.81*t,y(:,3)*10^7,12.81 * t,((y(:,4) + y(:,7))./((y(:,3) + y(:,6)))));
xlabel('Axial position,12.81*t, (m)'); 
ylabel(yax(1),'y_9,^{\circ}C, y3,mol/liter');%left y axis 
ylabel(yax(2),'\mu_n');%right y axis
% %plot 2  - Monomer conversion vs time)
% plot(12.81 * t,(1 - y(1:end,2)./(y(2,1))));
% xlabel('Axial position (m)');
% ylabel('Monomer conversion');
% figure;
% %plot 3 - Initiator conversion vs time)
% plot(12.81 * t,(y(1:end,1)./(y(1,1))));
% xlabel('Axial position (m)');
% ylabel('Intiator conversion');
% figure;
% %plot 4 - Free radical conversion vs time)
% plot(12.81 * t,y(1:end,3));
% xlabel('Axial position (m)');
% ylabel('Free radical conversion (mol/liter)');
% figure;
% %plot 5  - Number avg. chain length vs time)
% plot(12.81 * t,((y(1:end,4) + y(1:end,7))./((y(1:end,3) + y(1:end,6)))));
% xlabel('Axial position (m)');
% ylabel('Number avg. chain length');
% 
% 
% 
% 

