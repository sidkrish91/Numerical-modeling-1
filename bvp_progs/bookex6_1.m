% to solve example 6.1 from CRE with MATLAB
% lnK = lnA - E/R*(1/T) is the relation
% a plot of lnK vs 1/T yields slope = -E/R and intercept lnA
T=[313.0;319.0;323.0;328.0;333.0]; % temperatures
k=[0.00043;0.00103;0.00180;0.00355;0.00717]; % k value in (sec)^-1
lnk=log(k);
invT=[1./T];
x=polyfit(invT,lnk,1); % evaluates the co-efficients of polynomial.
A=exp(c(2);
E=-c(1)*8.314;

