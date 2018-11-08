function dy = AtmoHW2_4_solve(t,y )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
k1 = 9.971*10^-3;
k2 = 6.0*10^-34;
k3 = 1.9*10^-14;
M = 2.446*10^19;
O2 = 0.21*M;

dy(1,1) = k3*y(3)*y(4);
dy(2,1) = -k2*y(2)*O2*M;
dy(3,1) = -k3*y(3)*y(4);
dy(4,1) = k2*y(2)*O2*M-k3*y(3)*y(4);
% dy(1,1) = -k1*y(1)+k3*y(3)*y(4);
% dy(2,1) = k1*y(1)-k2*y(2)*O2*M;
% dy(3,1) = k1*y(1)-k3*y(3)*y(4);
% dy(4,1) = k2*y(2)*O2*M-k3*y(3)*y(4);
end

