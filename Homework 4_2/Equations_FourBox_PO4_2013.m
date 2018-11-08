function dcdt = Equations_FourBox_PO4_2013(t,c)

Definitions_ThreeBox_PO4_AtmCO2_2012;

%	Create the the 'A' Matrix
Apo4 = zeros(4,4);
Apo4(1,:) = [-fLH-fLD-T-kL*VL fLH fLD+T 0];
Apo4(2,:) = [fLH+T -fLH-fHD-T-kH*VH fHD 0];
Apo4(3,:) = [fLD+kL*VL fHD+T+kH*VH -fHD-fLD-T 0];
Apo4(4,:) = [0 0 0 0];
bpo4 = [0 0 0 0];

%   Evaluate the gradients for po4 in each box

next = Apo4*c(po4_1:po4_4).*InvOceanVolArray';
dcdt(po4_1) = sum(next(1)+bpo4(1));
dcdt(po4_2) = sum(next(2)+bpo4(2));
dcdt(po4_3) = sum(next(3)+bpo4(3));
dcdt(po4_4) = sum(next(4)+bpo4(4));
dcdt = dcdt';
