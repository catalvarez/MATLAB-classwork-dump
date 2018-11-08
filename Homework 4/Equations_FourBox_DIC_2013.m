function dDdt = Equations_FourBox_DIC_2013(t,D)

Definitions_ThreeBox_PO4_AtmCO2_2012;


%	Create the the 'A' Matrix

%DIC
ADIC = zeros(4,4);
ADIC(1,:) = [-fLH-fLD-T-kL*VL*rdic_po4 fLH fLD+T 0];
ADIC(2,:) = [fLH+T -fLH-fHD-T-kH*VH*rdic_po4 fHD 0];
ADIC(3,:) = [fLD+kL*VL*rdic_po4 fHD+T+kH*VH*rdic_po4 -fHD-fLD-T 0];
ADIC(4,:) = [0 0 0 0];
bDIC = [0 0 0 0];


%   Evaluate the gradients for DIC in each box

%DIC
next = ADIC*D(1:4).*InvOceanVolArray';
dDdt(1) = sum(next(1)+bDIC(1));
dDdt(2) = sum(next(2)+bDIC(2));
dDdt(3) = sum(next(3)+bDIC(3));
dDdt(4) = sum(next(4)+bDIC(4));
dDdt = dDdt';