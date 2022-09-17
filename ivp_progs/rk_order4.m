% to write a program to calculate ODE using range kutta 4th order technique
%function is dy/dx=(x+y)*sin(xy)
function y = rk_order4(h,guess) % please input step size and initial guess while executing the function
clc;
 %  h = step size
x = [0:h:2]; %calculates f_y upto x =2 in steps of h. please change the limit if needed 
y = zeros(1, length(x)); %y values are saved in a row vector with no. of columns as length(x) ie. one per step
y(1) = guess; % initial value
f_y = @(x,y) (x + y) * sin(x * y); % function declaration. function is declared for 2 variables. can change the function if needed

for i=1:(length(x)-1) %length(x)-1 due to initial guess already present for h=0 as guess = y(1)
    k_1 = f_y (x(i), y(i));
    k_2 = f_y (x(i) + h / 2, y(i) + k_1 * h / 2);
    k_3 = f_y (x(i) + h / 2, y(i) + k_2 * h / 2);
    k_4 = f_y (x(i) + h / 2, y(i) + k_3 * h / 2);
    y(i+1) = y(i) + h * ((k_1 * 1/ 6) + (k_2 * 1 / 3) + (k_3 * 1 / 3) + (k_4 * 1 / 6));
end
plot (x,y);

