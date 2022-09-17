%solution file for SKGupta_5_27_batch
clc
clear
%declaring glabal variables
global Ai_o; 
global Ei_o;
global Ai_c;
global delS;
global delH;
global Ei_c;

%initilizing global variables
 Ai_o = zeros(3,1);
 Ei_o = zeros(3,1);
 Ai_c = zeros(3,1);
 delS = zeros(3,1);
 delH = zeros(3,1);
 Ei_c = zeros(3,1);
 
 %retrieving values from excel file
 excel_data = xlsread('SKGupta_5.27-data.xlsx','A2:F4');
 
 %to parse values to global variables
Ai_o = excel_data(1:end,1);
Ei_o = excel_data(1:end,2);
Ai_c = excel_data(1:end,3);
Ei_c = excel_data(1:end,4);
delH = excel_data(1:end,5);
delS = excel_data(1:end,6);

%intializing input variables
y0 = [8.8, 0, 0, 0, 0.177];
tspan = [0 12];

%running ODE45
[t,y] = ode45(@SKGupta_5_27_batch, tspan, y0);

%plotting of curves

%plot 1 - caprolactin conversion vs time)
plot(t,(1-y(1:end,1)/y(1,1)));
xlabel('time (hr)');
ylabel('Caprolactan conversion');
figure;
%plot  - Average molecular weight vs time)
plot(t,113*y(1:end,4)./(y(1:end,3)));
xlabel(' time (hr)');
ylabel('Average molecular weight');


