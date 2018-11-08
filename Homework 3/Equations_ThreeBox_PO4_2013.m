function dcdt = Equations_ThreeBox_PO4_2013(t,c)

Definitions_ThreeBox_PO4_2013;

%	Create the the 'A' Matrix

Apo4 = zeros(3,3);
Apo4(1,:) = [-fLH-fLD-T-kL*VL fLH fLD+T];
Apo4(2,:) = [fLH+T -fLH-fHD-T fHD];
Apo4(3,:) = [fLD+kL*VL fHD+T -fHD-fLD-T];
% Apo4(1,:) = [-fLH-fLD-T-kL*VL fLH fLD+T];
% Apo4(2,:) = [fLH+T -fLH-fHD-T-kH*VH fHD];
% Apo4(3,:) = [fLD+kL*VL fHD+T+kH*VH -fHD-fLD-T];

bpo4 = [0 -1e-14 1e-14];

%   Evaluate the gradients for po4 in each box

next = Apo4*c(po4_1:po4_3).*InvOceanVolArray';

dcdt(po4_1) = sum(next(1)+bpo4(1));
dcdt(po4_2) = sum(next(2)+bpo4(2));
dcdt(po4_3) = sum(next(3)+bpo4(3));

testt = 2*t;

dcdt = dcdt';