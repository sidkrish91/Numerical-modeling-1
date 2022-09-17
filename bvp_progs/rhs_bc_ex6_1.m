function bc = rhs_bc_ex6_1( yL,yR )
%boundary conditions as f(x) = 0
Pe = 6;
bc = [yL(2) - Pe * (yL(1) - 1); yR(2)];


end

