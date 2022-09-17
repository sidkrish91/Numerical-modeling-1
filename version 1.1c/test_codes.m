%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global prop_data2;
Temp = 493.2
del_a1 = 1*prop_data2(1,9)+1*prop_data2(2,9)-1*prop_data2(4,9)-3*prop_data2(5,9); 
del_a2 = 1*prop_data2(3,9)+1*prop_data2(2,9)-1*prop_data2(5,9)-1*prop_data2(4,9);
del_b1 = 1*prop_data2(1,10)+1*prop_data2(2,10)-1*prop_data2(4,10)-3*prop_data2(5,10);
del_b2 = 1*prop_data2(3,10)+1*prop_data2(2,10)-1*prop_data2(5,10)-1*prop_data2(4,10);
del_c1 = 1*prop_data2(1,11)+1*prop_data2(2,11)-1*prop_data2(4,11)-3*prop_data2(5,11); 
del_c2 = 1*prop_data2(3,11)+1*prop_data2(2,11)-1*prop_data2(5,11)-1*prop_data2(4,11);
del_d1 = 1*prop_data2(1,12)+1*prop_data2(2,12)-1*prop_data2(4,12)-3*prop_data2(5,12); 
del_d2 = 1*prop_data2(3,12)+1*prop_data2(2,12)-1*prop_data2(5,12)-1*prop_data2(4,12);
del_e1 = 1*prop_data2(1,13)+1*prop_data2(2,13)-1*prop_data2(4,13)-3*prop_data2(5,13); 
del_e2 = 1*prop_data2(3,13)+1*prop_data2(2,13)-1*prop_data2(5,13)-1*prop_data2(4,13);

T1 = 298;
fun1 = @(T) del_a1+del_b1*T+del_c1*T.^2+del_d1*T.^3+del_e1*T.^4;
fun2 = @(T) del_a2+del_b2*T+del_c2*T.^2+del_d2*T.^3+del_e2*T.^4;

cp1 = integral(fun1,T1,Temp)*8.314;%j/mol/K
cp2 = integral(fun2,T1,Temp)*8.314;%j/mol/K
del_h1 = -49500+cp1 %MeOH 
del_h2 = 41000+cp2; %RWGS

delH1 = -49500 + (del_a1*(Temp-T1) + del_b1*(Temp^2-T1^2)/2 + del_c1*(Temp^3-T1^3)/3 + ...
    del_d1*(Temp^4-T1^4)/4 + del_e1*(Temp^5-T1^5)/5)*8.314