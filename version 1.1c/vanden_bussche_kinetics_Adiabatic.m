function dy_dt = vanden_bussche_kinetics_Adiabatic(~,y)
% To solve a system of coupled ODE's for the mass balance for a
%Fixed bed reactor. Adiabatic reactor conditions were modeled. Unsteady
%state model coded

dy_dt=zeros(6,1);
%%%%%%%%%%%
global kin_data;
global p_inlet;
global u_gas;
global prop_data2;
global rho_gas_mix;
global t_inlet;

rho_bed = 1775; %density of catalyst bed, kg-cat/m^3
%%%%%%%%
p_par=zeros(5,1);
% T=y(7);
c_tot = p_inlet*10^5/8.314/t_inlet; %concentration of gas mixture, mol/m^3
T = y(6);
%defining the rate constants
k=zeros();

for j=1:5 %5 reaction rate constants
    A = kin_data(j,1);
    B = kin_data(j,2);
    k(j) = A * exp(B/8.314/T); %rate constant 8.314J/kelvin/mole,  
end

%equilibrium rate constants
keq1= 10^(3066/T - 10.592);
keq2 = (10^(-2073/T + 2.029))^-1;

for i=1:5
    p_par(i) = p_inlet*y(i);
end
    
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
% if r1<r2
%     disp('rwgs is higher')
% else
%     disp('MeoH rate is higher')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dy_dt(1) = 0.5*rho_bed * (r1) / (c_tot*u_gas); %change of methanol mole fraction with axial length                 %
dy_dt(2) = 0.5*rho_bed * (r1+r2)/ (c_tot*u_gas); %change of water mole fraction with axial length                 %
dy_dt(3) = 0.5*rho_bed * (r2)/ (c_tot*u_gas); %change of CO mole fraction with axial length                       %
dy_dt(4) = 0.5*rho_bed * (-r1-r2)/ (c_tot*u_gas); %change of CO2 mole fraction with axial length                  %
dy_dt(5) = 0.5*rho_bed * (-3*r1-r2)/ (c_tot*u_gas); %change of H2 mole fraction with axial length                 %
% dy_dt(6) = 0;
% 1-sum(dy_dt(5,:)); %change of N2 mole fraction with axial length                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp_pure2 = zeros();
for i=1:5
    cp_pure2(i) = 8.314*(prop_data2(i,9) + prop_data2(i,10)*T + prop_data2(i,11)*T^2 +...
        prop_data2(i,12)*T^3 + prop_data2(i,13)*T^4)*10^3/prop_data2(i,1);
    % NASA Coefficients
end

%specific heat mixture J/kg*kelvin 
cp_mix = 0; 
for i=1:length(cp_pure2)
cp_mix = cp_mix + y(i)*cp_pure2(i); %specific heat of mixture based on mole fraction 
end

del_h = enthalpy_calc(T);
% heat balance
dy_dt(6) = rho_bed*(-r1*del_h(1)-r2*del_h(2))/(rho_gas_mix*cp_mix*u_gas); 
 
end

