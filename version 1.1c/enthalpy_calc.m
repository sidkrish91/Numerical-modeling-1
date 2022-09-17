function del_h = enthalpy_calc( Temp )
%FOR VERSION 1 CODE
global prop_data2;
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

T1 = 298; %reference temperature, 298K
fun1 = @(T) del_a1+del_b1.*T+del_c1.*T.^2+del_d1.*T.^3+del_e1.*T.^4;
fun2 = @(T) del_a2+del_b2.*T+del_c2.*T.^2+del_d2.*T.^3+del_e2.*T.^4;

cp1 = integral(fun1,T1,Temp)*8.314;%j/mol/K * K
cp2 = integral(fun2,T1,Temp)*8.314;%j/mol/K * K
del_h(1) = -90.5e3+cp1; %MeOH 
del_h(2) = 41.0e3+cp2; %RWGS
end

