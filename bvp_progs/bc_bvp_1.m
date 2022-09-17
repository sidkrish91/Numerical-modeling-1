function bc = bc_bvp_1( yL,yR)
%pass the extreme boundary conditions or intermediate BC's
bc = [yL(1)-3;yR(1)+2*yR(2)-5];

end
