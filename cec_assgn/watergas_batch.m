%function program to solve system of ODE's for watergas reaction from
%ethanol
function dy=watergas_batch(t,y)

%parsing intial values (molar conc.)
ycd=y(1);  %carbon dioxide
ycm=y(2);  %carbon monoxide
yh=y(3);   % hydrogen
ym=y(4);   % methane
yw=y(5);   % water
ye=y(6);   %ethanol

%intializing variables
dy=zeros(6,1);
ptot=1; %total pressure

%reaction rate constants
k1=2.9e-5;
k2=3.1e-4;
k3=1e-4;
k4=2.4e-4;
Kw=37.4;
Ke=61.7;
Km=1135;
K3=0.0842;
K4=0.01434;

%partial prressures
Pe=ptot*ye;
Pw=ptot*yw;
Pm=ptot*ym;
Pcm=ptot*ycm;
Ph=ptot*yh;
Pcd=ptot*ycd;

%defining rate equations
DEN=(1+Pe*Ke+Pw*Kw+Pm*Km);
r1=k1*Ke*Pe/DEN*10^3;
r2=(k2*Ke*Kw*Pe*Pw)/(DEN^2)*10^3;
r3=k3*Km*Kw*(Pm*Pw-(Pcm*Ph^3/K3))/DEN^2*10^3;
r4=k4*Km*Kw*(Kw*Pm*Pw^2-(Pcd*Ph^4/K4))/DEN^3*10^3;

%govering differential equations
dy(1)=r2+r4;
dy(2)=r1+r3;
dy(3)=r1+2*r2+3*r3+4*r4;
dy(4)=r1+r2-r3-r4;
dy(5)=-1*r3-2*r4-r2;
dy(6)=-r1-r2;
