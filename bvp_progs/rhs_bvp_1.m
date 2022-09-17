function rhs = rhs_bvp_1(x,y)
%contains the system of first order differntial equations. The parent 2nd
%order equation was broken down into sytem of 1st order differential
%equations 

rhs = [y(2); 5-3*y(2)-6*y(1)]; 

end


