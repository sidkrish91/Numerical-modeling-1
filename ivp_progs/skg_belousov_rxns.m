function dy_dt = skg_belousov_rxns(~,y)
%this contains the time dependent mass balances of Belousov reactions.
%a mixture of 1.25M H2SO4, 0.0125m KBrO3, 0.001m Ce(NH4)2(NO3)2, 0.025M
%CH2(COOH)2
dy_dt=zeros(3,1);
dy_dt = [77.27*(y(2)-y(1)*y(2)+y(1)-pow2(y(1))*8.375e-6); (-y(2)-y(1)*y(2)+y(3))/77.27; 0.161*(y(1)-y(3))];
end

