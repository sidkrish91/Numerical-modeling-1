function  dy_dt = chemburkar_func(t,y,gamma)
%to solve the Chemburkar et. al which deals with dynamics of a CSTR
% reference: SK Gupta page 201

dy_dt = zeros(3,1);
% defining parameters



dy_dt(1) = 1 - y(1) - gamma(1)*y(1)*exp(gamma(2)*y(3)/(gamma(2) + y(3)));
dy_dt(2) = -1*y(2) + gamma(1)*y(1)*exp(gamma(2)*y(3)/(gamma(2) + y(3))) - ...
    gamma(1)*gamma(3)*y(2)*exp(gamma(4)*gamma(2)*y(3)/(gamma(2) + y(3)));
dy_dt(3) = -1*y(3) - gamma(5)*(y(3) - gamma(6)) + gamma(1)*gamma(7)*y(1)*exp(gamma(2)*y(3)/(gamma(2) + y(3)))+...
    gamma(1)*gamma(7)*gamma(3)*gamma(8)*y(2)*exp(gamma(4)*gamma(2)*y(3)/(gamma(2) + y(3)));
end

