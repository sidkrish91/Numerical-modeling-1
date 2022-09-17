function dy_dz = vanden_bussche_kinetics_fbrode(~,y)
% To solve a system of coupled ODE's for the mass and energy balance for a
%Fixed bed reactor. NINA reactor equations used 

dy_dz=zeros(7,1);

%%%%%%%%%%%
global kin_data;
global p_inlet;
global u_gas;
global id_reactor;
global rho_gas_mix;
global cp_mix;
%%%%%%%%%%%

rho_bed = 1140; %density of catalyst bed, kg-cat/m^3

% temperature variable
T = y(7);

% Partial pressure of components
p_par = zeros();
% note, p_par order: MeOH,H2O,CO,CO2,H2,N2-6 components
for j = 1:6
    p_par(j) = p_inlet * y(j); %order: Methanol, water, CO, CO2, H2, N2
end

c_tot = p_inlet*10^5/8.314/493.2; %concentration of gas mixture, mol/m^3

%defining the rate constants
k=zeros();

for j=1:5 %5 reaction rate constants
    A = kin_data(j,1);
    B = kin_data(j,2);
    k(j) = A * exp(B/(8.314*T)); %rate constant 8.314J/kelvin/mole,  
end

%equilibrium rate constants
keq1= 10^(3066/T - 10.592);
keq2 = (10^(-2073/T + 2.029))^-1;

%rate equations
den = 1 + k(2)*(p_par(2)/p_par(5)) + k(3)*(p_par(5)^0.5) + k(4)*p_par(2);

a = k(1)*p_par(4)*p_par(5);
b =   1 - p_par(2)*p_par(1)/(keq1*p_par(4)*(p_par(5)^3));
r1 = a*b/den^3;%methanol formation rate constant
% mol/kg.cat*sec

a = k(5)*p_par(4);
b = 1-keq2*p_par(2)*p_par(3)/(p_par(4) * p_par(5));
r2 = a*b/den;
%reverse water gas shift reaction rate constant  mol/kg.cat*sec
if r1<r2
    disp('rwgs is higher')
else
    disp('MeoH rate is higher')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy_dz(1) = rho_bed * (r1) / (c_tot*u_gas); %change of methanol mole fraction with axial length                 %
dy_dz(2) = rho_bed * (r1+r2)/ (c_tot*u_gas); %change of water mole fraction with axial length                 %
dy_dz(3) = rho_bed * (r2)/ (c_tot*u_gas); %change of CO mole fraction with axial length                       %
dy_dz(4) = rho_bed * (-r1-r2)/ (c_tot*u_gas); %change of CO2 mole fraction with axial length                  %
dy_dz(5) = rho_bed * (-3*r1-r2)/ (c_tot*u_gas); %change of H2 mole fraction with axial length                 %
dy_dz(6) = 0; %change of N2 mole fraction with axial length                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% heat balance
dy_dz(7) = (rho_bed*(r1*49430-r2*41000)-4*18000*(T-550)/id_reactor)/(rho_gas_mix*cp_mix); 
 
end

