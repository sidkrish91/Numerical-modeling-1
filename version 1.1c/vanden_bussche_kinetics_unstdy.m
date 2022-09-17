function dy_dt = vanden_bussche_kinetics_unstdy(~,y)
% To solve a system of coupled ODE's for the mass balance for a
%Fixed bed reactor. Adiabatic reactor conditions were modeled. Unsteady
%state model coded

dy_dt=zeros(600,1);

%%%%%%%%%%%
global kin_data;
global p_inlet N;
global u_gas dz;
global prop_data2;
global rho_gas_mix;
global t_inlet;

rho_bed = 1775; %density of catalyst bed, kg-cat/m^3
%%%%%%%%%%%  
for k=1:6 %compartmentalising the nodes
    switch k
        case 1 %for meoh
            y_meoh=y(k:N*k); %mole fraction of component
            p_meoh=y(k:N*k).*p_inlet; %partial pressure of component
        case 2 %for H2O
            y_water=y((k-1)*N+1:N*k);
            p_water=y((k-1)*N+1:N*k).*p_inlet;
        case 3 %for CO
            y_CO=y((k-1)*N+1:N*k);
            p_CO=y((k-1)*N+1:N*k).*p_inlet;
        case 4 %for CO2
            y_CO2=y((k-1)*N+1:N*k);
            p_CO2=y((k-1)*N+1:N*k).*p_inlet;
        case 5 %for H2
            y_H2=y((k-1)*N+1:N*k);
            p_H2=y((k-1)*N+1:N*k).*p_inlet;
%         case 6 %for N2
%             y_N2=y((k-1)*N+1:N*k);
%             p_N2=y((k-1)*N+1:N*k).*p_inlet;
        case 6 % for temperature
            T=y((k-1)*N+1:N*k);
    end
end

%%. sectional code to saying change wrt time at inlet and exit of node is zero.
% for i = 1:N %iterating with the nodes
%     for m=1:length(dy_dt)
%         if rem(i,N)==1 %at inlet node, no gradient with respect to space
%             dy_dt(i) = 0; 
%     elseif rem(i,N)==0 %for all end nodes
%             dy_dt(i) = 0;
%         end
%     end
% end

c_tot = p_inlet*10^5/8.314/t_inlet; %concentration of gas mixture, mol/m^3
%% all the FDM calculations begin from the 1st array element, until N node
 
%defining the rate constants

k=zeros(1,5);
 %5 reaction rate constants
    for m = 1:5
        A = kin_data(m,1);
        B = kin_data(m,2);
        k(m) = A * exp(B/8.314/T);%rate constant 8.314J/kelvin/mole, 
    end
%equilibrium rate constants
keq1= 10^(3066/T - 10.592);
keq2 = (10^(-2073/T + 2.029))^-1;

for i=1:N %iterations with Nodes as multiples of 10. For simplicity

%defining the rate constants

k=zeros(1,5);
 %5 reaction rate constants
    for m = 1:5
        A = kin_data(m,1);
        B = kin_data(m,2);
        k(m) = A * exp(B/8.314/T(i));%rate constant 8.314J/kelvin/mole, 
    end
%equilibrium rate constants
keq1= 10^(3066/T(i) - 10.592);
keq2 = (10^(-2073/T(i) + 2.029))^-1;

    
    %rate equations
        den = 1 + k(2)*(p_water(i))/(p_H2(i)) + k(3)*((p_H2(i))^0.5)...
    + k(4)*(p_water(i));

a = k(1)*(p_CO2(i))*(p_H2(i));
b =   1 - (p_water(i))*(p_meoh(i))/...
    (keq1*(p_CO2(i))*((p_H2(i))^3));
r1 = a*b/den^3;%methanol formation rate constant
% mol/kg.cat*sec

a = k(5)*(p_CO2(i));
b = 1-keq2*(p_water(i))*(p_CO(i))/...
    ((p_CO2(i)) * (p_H2(i)));
r2 = a*b/den;
%reverse water gas shift reaction rate constant  mol/kg.cat*sec
% if r1<r2
%     disp('rwgs is higher')
% else
%     disp('MeoH rate is higher')
% end
%points to note: 1) at inlet and end nodes the dy_dt=0. hence inlet & exit
%nodes of all species is zero     
if i==1
for j = 1:6   %stepping wrt species & temp           
        dy_dt((j-1)*N+i)=0; %change wrt time at inlet node =0
end
elseif rem(i,N)==0  
    for j = 1:6   %stepping wrt species & temp        
        dy_dt((j-1)*N+i)=0; %change wrt time at inlet node =0
    end
else
    for j = 1:6 
        switch j
            case 1 %methanol mass balance at ith node
               dy_dt((j-1)*N+i) = -u_gas*(y_meoh(i)-y_meoh(i-1))/dz + 0.5*rho_bed * (r1) / (c_tot); 
                %change of methanol mole fraction with axial length                 
                
        case 2 %water fraction  
            dy_dt((j-1)*N+i) = -u_gas*(y_water(i)-y_water(i-1))/dz + 0.5*rho_bed * (r1+r2) / (c_tot);
                 %change of water mole fraction with axial length                 %
            case 3 %CO fraction
                dy_dt((j-1)*N+i)=-u_gas*(y_CO(i)-y_CO(i-1))/dz + 0.5*rho_bed * (r2) / (c_tot);
                %change of CO mole fraction with axial length                       
         
        case 4 %CO2 fraction
            dy_dt((j-1)*N+i)=-u_gas*(y_CO2(i)-y_CO2(i-1))/dz + 0.5*rho_bed * (-r1-r2) / (c_tot);
                %change of CO2 mole fraction with axial length
                
            
            case 5 %hydrogen fraction
                dy_dt((j-1)*N+i)=-u_gas*(y_H2(i)-y_H2(i-1))/dz + 0.5*rho_bed * (-3*r1-r2) / (c_tot);
                %change of H2 mole fraction with axial length
            
%             case 6
%                 dy_dt((j-1)*N+i) = 0; %change of N2 mole fraction with axial length
            case 6
                dy_dt((j-1)*N+i)
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % cp_pure2 = zeros();
% % % for i=1:6
% % %     cp_pure2(i) = 8.314*(prop_data2(i,9) + prop_data2(i,10)*T + prop_data2(i,11)*T^2 +...
% % %         prop_data2(i,12)*T^3 + prop_data2(i,13)*T^4)*10^3/prop_data2(i,1);
% % %     % NASA Coefficients
% % % end
% % % 
% % % %specific heat mixture J/kg*kelvin 
% % % cp_mix = 0; 
% % % for i=1:length(cp_pure2)
% % % cp_mix = cp_mix + y(i)*cp_pure2(i); %specific heat of mixture based on mole fraction 
% % % end
% % % 
% % % del_h = enthalpy_calc(T);
% % % % heat balance
% % % dy_dt(7) = rho_bed*(-r1*del_h(1)-r2*del_h(2))/(rho_gas_mix*cp_mix*u_gas); 
 
end
end
end

