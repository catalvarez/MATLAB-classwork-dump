function dy = HW3_P2_solve (t,y)
jO3O1D = 2.7*10^-4;
kO1DCH4 = 1.5*10^-10;
kO1DH2O = 2.2*10^-10;
kO1DO2 = 3.2*10^-11*exp(70/227);
kO1DN2 = 1.8*10^-11*exp(110/227);
H2O = 1.68*10^12;
CH4 = 4.2*10^11;
O3 = 4.998*10^12;
M = 4.2*10^17;
N2 = M*.79;
O2 = M*.21;

dy(1) = jO3O1D*O3 - kO1DCH4*y(1)*CH4 - kO1DH2O*y(1)*H2O - kO1DO2*y(1)*O2 - kO1DN2*y(1)*N2; %O(1D)

end