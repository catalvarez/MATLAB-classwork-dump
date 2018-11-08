function dAdt = Equations_FourBox_Alk_2013(t,A)

Definitions_FourBox_PO4_AtmCO2_2013;

%	Create the the 'A' Matrix
%Alk
AAlk = zeros(4,4);
AAlk(1,:) = [-fLH-fLD-T-kL*VL*ralk_po4 fLH fLD+T 0];
AAlk(2,:) = [fLH+T -fLH-fHD-T-kH*VH*ralk_po4 fHD 0];
AAlk(3,:) = [fLD+kL*VL*ralk_po4 fHD+T+kH*VH*ralk_po4 -fHD-fLD-T 0];
AAlk(4,:) = [0 0 0 0];
bAlk = [0 0 0 0];

%   Evaluate the gradients in each box

%Alk
next = AAlk*A(1:4).*InvOceanVolArray';
dAdt(1) = sum(next(1)+bAlk(1));
dAdt(2) = sum(next(2)+bAlk(2));
dAdt(3) = sum(next(3)+bAlk(3));
dAdt(4) = sum(next(4)+bAlk(4));
dAdt = dAdt';