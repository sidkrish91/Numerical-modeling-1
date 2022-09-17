% a program depicting the polymerization of nylon 6 in a well mixed
% isothermal batch reactor (no vaporization)

function dy = SKGupta_5_27_batch(t,y)

%declaring variables (global variables for ease of use)

global Ai_o; 
global Ei_o;
global Ai_c;
global delS;
global delH;
global Ei_c;

%intializing variables & declaring constants
dy=zeros(5,1);
k = zeros(3,1); %reaction rate constant
ko = zeros(3,1); %reaction rate constant
T = 538;

%parsing intial values to respective variables
y1 = y(1);
y2 = y(2);
y3 = y(3);
y4 = y(4);
y5 = y(5);


%passing values to reaction rate constants using global variables
for i = 1:length(k)
    k(i) = Ai_o(i) * exp(( -1 * Ei_o(i))/(1.987 * T)) + Ai_c(i) * y3 * exp((-1 * Ei_c(i))/(1.987 * T));
    ko(i) = k(i) * inv(exp((delS(i)/1.987) - (delH(i)/(1.987 * T))));
end

%polmerization reaction rate equations
dy(1) = -1 * (k(1) * y1 * y5) + (ko(1) * y2) - (ko(3) * (y3 - y2));
dy(2) = (k(1) * y1 * y5) - (ko(1) * y2) - (2 * k(2) * y2 * y3) + (2 * ko(2) * y5 * (y3 - y2)) - (k(3) * y1 * y2) + (ko(3) * y2);
dy(3) = (k(1) * y1 * y5) - (ko(1) * y2) - (k(2) * (y(3)^2)) + (ko(2) * y5 * (y4 - y3));
dy(4) = (k(1) * y1 * y5) - (ko(1) * y2) + (k(3) * y3 * y1) - (ko(3) * (y3 - y2));
dy(5) = (-1 * k(1) * y1 * y5) + (ko(1) * (y2)) + (k(2) * (y3^2))- (ko(2) * y5 * (y4 - y5));
