% a program depicting the high pressure polymerization of ethylene in 
% plug flow reactor 

function dy = SKGupta_5_12_plug(t,y)

%intializing variables
dy=zeros(9,1);

%parsing intial values to respective variables
y1 = y(1);
y2 = y(2);
y3 = y(3);
y4 = y(4);
y5 = y(5);
y6 = y(6);
y7 = y(7);
y8 = y(8);
y9 = y(9);

%declaring constants and reaction rates
kL = (1.6*10^16)*exp(-38400/(1.987*(y9+273)));
kP = (5.888*10^7)*exp(-5987/(1.987*(y9+273)));
kT = (1.3963*10^9)*exp(354/(1.987*(y9+273)));
beta1 = 76.9231;
beta2 = 0.06651;
beta3 = 180;

% defining reaction rate equation

dy(1) = -kL*y1;
dy(2) = -kP*y2*y3;
dy(3) = (2*kL*y1) - (kT*(y3^2));
dy(4) = (2*kL*y1) + (kP*y2*y3) - (kT*y3*y4);
dy(5) = (2*kL*y1) + (kP*y2*(2*y4+y3)) - (kT*y3*y5);
dy(6) = (0.5*kT*(y3^2));
dy(7) = (kT*y3*y4);
dy(8) = kT*(y3*y5+y4^2);
dy(9) = (beta1*kP*y2*y3) + beta2*(beta3-y9);







