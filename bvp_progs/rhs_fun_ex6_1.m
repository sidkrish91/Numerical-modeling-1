function rhs = rhs_fun_ex6_1( x,y )
%the governing equations in terms of 1st order differential equations
Pe = 6; %peclet no
Da = 2; %damkohler no
rhs = [y(2); (Da*(y(1)^2) + y(2))*Pe];
end

