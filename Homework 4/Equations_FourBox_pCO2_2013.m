function dpdt = Equations_FourBox_pCO2_2013(t,p)

Definitions_ThreeBox_PO4_AtmCO2_2012;

%	Create the the 'A' Matrix

%DIC
ApCO2 = zeros(4,4);
ApCO2(1,:) = [-fLH-fLD-T-kL*VL*rdic_po4-PV_1*Area_1*pco2 fLH fLD+T PV_1*Area_1*B_1*pco2];
ApCO2(2,:) = [fLH+T -fLH-fHD-T-kH*VH*rdic_po4-PV_2*Area_1*pco2 fHD PV_2*Area_2*B_2*pco2];
ApCO2(3,:) = [fLD+kL*VL*rdic_po4 fHD+T+kH*VH*rdic_po4 -fHD-fLD-T 0];
ApCO2(4,:) = [PV_1*B_2*Area_1*pco2 PV_2*B_2*Area_2*pco2 0 PV_1*Area_1*pco2+PV_2*Area_2*pco2];
bpCO2 = [0 0 0 0];

%   Evaluate the gradients for pCO2 in each box

%pCO2
next = ApCO2*p(1:4).*InvOceanVolArray';
dpdt(1) = sum(next(1)+bpCO2(1));
dpdt(2) = sum(next(2)+bpCO2(2));
dpdt(3) = sum(next(3)+bpCO2(3));
dpdt(4) = sum(next(4)+bpCO2(4));
dpdt = dpdt';