function dy = HW3_NOx_solve (t,y)

kO2 = 6.0*10^-11;
kO2OM = 6.0*10^-34;
jO3 = 0.0012;
kO3O = 9.16*10^-16;
kHO2O3 = 1.0e-14*exp(-490/227);
kHO2O = 3.0e-11*exp(200/227);
kHO2OH = 4.8e-11*exp(250/227);
jO3O1D = 2.7*10^-4;
kO1DCH4 = 1.5*10^-10;
kO1DH2O = 2.2*10^-10;
kO1DO2 = 3.2*10^-11*exp(70/227);
kO1DN2 = 1.8*10^-11*exp(110/227);
kNO2OH = 3.1915e-13;
jHNO3 = 1.5e-5;
kOHO3 = 1.7*10^-12*exp(-940/227);
kHO2NO = 1.05e-11;
kO3NO = 4.05e-15;
kO3NO2 = 2.465e-18;
jNO2 = 1.5e-2;
kN2OO1D = 6.7e-11;
kN2OO1D2 = 4.9e-11;
kNO2O = 1.24e-11;
jN2O = 5e-8;

H2O = 1.68*10^12;
CH4 = 4.2*10^11;
M = 4.2*10^17;
N2 = M*.79;
O2 = M*.21;
HNO3 = 9.4643e8;
N2O = 155e-9*M;


dy(1,1) = kO2OM*O2*y(2)*M -jO3*y(1) -kO3O*y(2)*y(1) -kOHO3*y(1)*y(3) -kHO2O3*y(1)*y(4) -kO3NO*y(6)*y(1) -kO3NO2*y(7)*y(1); %O3
dy(2,1) = 2*kO2*O2 -kO2OM*O2*y(2)*M +jO3*y(1) -kO3O*y(2)*y(1) -kHO2O*y(4)*y(2) +jNO2*y(7) -kNO2O*y(7)*y(2); %O
dy(3,1) = -kOHO3*y(3)*y(1) +kHO2O*y(4)*y(2) +kHO2O3*y(4)*y(1) -kHO2OH*y(3)*y(4) +kO1DCH4*y(5)*CH4 +2*kO1DH2O*y(5)*H2O -kNO2OH*y(7)*y(3) +jHNO3*HNO3 +kHO2NO*y(4)*y(6);  %OH
dy(4,1) = -kHO2OH*y(3)*y(4) +kOHO3*y(3)*y(1) -kHO2O3*y(4)*y(1) -kHO2O*y(4)*y(2) -kHO2NO*y(4)*y(6) ;   %HO2
dy(5,1) = jO3O1D*y(1) -kO1DCH4*y(5)*CH4 -kO1DH2O*y(5)*H2O -kO1DO2*y(5)*O2 -kO1DN2*y(5)*N2 -kN2OO1D*N2O*y(5) -kN2OO1D2*N2O*y(5) +jN2O*N2O;   %O1D
dy(6,1) = -kHO2NO*y(4)*y(6) -kO3NO*y(6)*y(1) +jNO2*y(7) +2*kN2OO1D*N2O*y(5) +kNO2O*y(7)*y(2);   %NO
dy(7,1) = -kNO2OH*y(7)*y(3) +jHNO3*HNO3 +kHO2NO*y(4)*y(6) +kO3NO*y(6)*y(1) -kO3NO2*y(7)*y(1) -kNO2O*y(7)*y(2) -jNO2*y(7);   %NO2

end