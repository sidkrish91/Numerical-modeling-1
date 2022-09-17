function rhs = rhs_bvp(x,y,beta)
%contains the system of first order differntial equations. The parent 2nd
%order equation was broken down into sytem of 1st order differential
%equations 

rhs = [y(2); (beta-100)*y(1)-10*y(1)^3]; %non linear differential equation

end

