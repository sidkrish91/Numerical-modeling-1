function bc = bc_bvp( yL,yR,beta )
%pass the extreme boundary conditions or intermediate BC's
bc = [yL(1);yL(2)-1;yR(1)];

end

