% This matlab file is used to run the function file to solve a
% pseudohomogeneous 1d model of a fixed bed reactor. VERSION 1.1c FOLDER.
% Unsteady state model for adiabatic condition
clc;clear;close all;
%defining parameters

global p_inlet;
global u_gas N;
global id_reactor;
global rho_gas_mix;
global kin_data;
global prop_data2;
global t_inlet y_inlet dz;
% global rate_MeoH;
% global rate_RWGS;
% rate_MeoH = 0;
% rate_RWGS = 0;

p_inlet = 50; %inlet pressure, Bar
t_inlet = 273.2+220; %inlet temperature, K
m_inlet = 2.8e-5; %inlet gas mass flow rate, kg/s
id_reactor = 0.016; %reactor diameter,m
L_reactor = 0.15; %reactor length, m

kin_data = xlsread('kinetics_data.xlsx','B2:C6'); 
% save('kinetic_data','kin_data'); %saving the excel data for use in matlab workspace
% Retrieving values from excel file
%  prop_data2 = xlsread('property data - with nasa coeffs.xlsx','D2:P7');
%  save('property_data2','prop_data2'); %saving the excel data 
load('kinetic_data.mat'); 
load('property_data2.mat');%Loading the property data with NASA coefficients
% load('property_data.mat')
%%%%%
y_inlet = [0.00 0.00 0.1 0.2 0.7]; 
%methanol, water, carbon monoxide, carbon dioxide, hydrogen and inert nitrogen
%%
%To calculate gas mixture molecular weight at inlet (g/mol)
mol_wt_mix = 0;
for i=1:length(y_inlet)
    mol_wt_mix = mol_wt_mix + y_inlet(i)* prop_data2(i,1);
end
%% density of gas mixture (kg/m^3)
rho_gas_mix = p_inlet * mol_wt_mix*10^2/ (8.314 * t_inlet);
%% to calculate superficial gas velocity
ac = 0.25*pi*id_reactor^2; %cross sectional area of reactor
u_gas = m_inlet/rho_gas_mix/ac;%velocity of inlet gas based on mass flow rate, m/s
%% To calculate the specific heat of pure gas component at the input temperature 
% and to return specific heat of individual gas component
%J/kg*kelvin, Order: Methanol,Water CO, CO2, H2 cp_pure = zeros(1,6); 
% for i=1:6
%     cp_pure(i) = prop_data(i,9) + prop_data(i,10)*t_inlet + prop_data(i,11)*t_inlet^2 +...
%         prop_data(i,12)*t_inlet^3 + prop_data(i,13)*t_inlet^4 +...
%         prop_data(i,14)*t_inlet^5 + prop_data(i,15)*t_inlet^6;
%     Specific heat based on polynomial correlations, Source:Calhoun,
%     Institutional Archive % of Naval Postgraduate School
% end
%
% cp_pure2 = zeros();
% for i=1:6
%     cp_pure2(i) = 8.314*(prop_data2(i,9) + prop_data2(i,10)*t_inlet + prop_data2(i,11)*t_inlet^2 +...
%         prop_data2(i,12)*t_inlet^3 + prop_data2(i,13)*t_inlet^4)*10^3/prop_data2(i,1);
%     % NASA Coefficients
% end
% %specific heat mixture J/kg*kelvin `  1
% cp_mix = 0; 
% for i=1:length(cp_pure2)
% cp_mix = cp_mix + y_inlet(i)*cp_pure2(i); %specific heat of mixture based on mole fraction 
% end

%%
z_span = [0 L_reactor];
y_span = [y_inlet t_inlet];
[z,y] = ode15s('vanden_bussche_kinetics_Adiabatic',z_span,y_span);
%% plotting the data
plot(z./L_reactor,y(1:end,6));
xlabel('Dimensionless length, z/L');
ylabel('Temperature, K');
figure;
plot(z./L_reactor,y(1:end,1:4).*100);
xlabel('Dimensionless length, z/L');
ylabel('Concentration, mole%');
legend('Methanol','Water','Carbon-monoxide','Carbon dioxide');
% figure;
% plot(z./L_reactor,y(1:end,5).*100);
% xlabel('Dimensionless length, z/L');
% ylabel('Concentration, mole%');
% legend('Hydrogen');

%% unsteady state 
N=100; %number of nodes, enter multiples of 10, preferably

% initializing the nodes at time t=0. Input nodes multiples of 10
y_initial=zeros(N*length(y_inlet),1);
%%initializing component mole fractions
% inserting 1e-4 as default values between node 1 and node 10
for i=1:length(y_initial)
        y_initial(i)=1e-3;
end
for m=1:length(y_inlet)
    y_initial(m*N) = y(end,m); %inserting steady state values at end nodes
    if m==1 %inserting inlet mole fractions of components
        y_initial(m)=y_inlet(m); 
    else
        y_initial((m-1)*N+1)=y_inlet(m);
    end
end
% %initializing temperatures at nodes
temp_initial=zeros(N,1);
% inserting 298k as default values between node 1 and node 10
for i=1:N
        temp_initial(i)=298;
end
temp_initial(1)=t_inlet;
temp_initial(N)=y(end,6);

% discretizing the nodes
dz=L_reactor/N;
t_span = [0 500]; %time stepping, seconds
y_span=cat(1,y_initial,temp_initial); %complete discretized grid, 1)methanol...
... 2)water, 3)CO, 4)CO2, 5)H2, 6)N2
% [t,yt]=ode15s('vanden_bussche_kinetics_unstdy',t_span,y_span);

%%
for k=1:5:length(t)
plot(yt(k,1:100));
ylim([0 0.1])
xlabel('reactor length')
ylabel('mole fraction')
text(5,0.09,sprintf('t=%1.2f',t(k)));
frame=getframe(1);
end