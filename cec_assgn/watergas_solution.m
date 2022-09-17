%solution code to run watergas_batch code 
clc
clear
tspan=[0:0.1:0.6]; %units consistency (for 10^-3 hrs)
y0=[0,0,0,0,0.09,0.016];
[t,y]=ode45(@watergas_batch,tspan,y0);
plot(t,y(1:end,1));
hold on
plot(t,y(1:end,2),'b--o');
hold on
plot(t,y(1:end,3),'c-*');
hold on
plot(t,y(1:end,4));
hold on
plot(t,y(1:end,5),'-.');
hold on
plot(t,y(1:end,6),'--');
xlabel('Space time(mgCat.min/(1000mol))');
ylabel('Molar function');
legend('y = Carbon Dioxide','y = Carbon monoxide','y = Hydrogen','y = Methane','y = Water','y = Ethanol','Location','northeast')