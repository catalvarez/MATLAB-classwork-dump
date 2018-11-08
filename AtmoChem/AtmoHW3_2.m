%% Problem 1
% [t,y] = ode15s(@HW3_O3_solve, [0 1e7], [0 0]);
% plot(t,y(:,1))
% O3SS = y(length(y),1)
% figure (2)
% plot(t,y(:,2))

%% Problem 2 a
% [t,y] = ode15s(@HW3_P2_solve, [0 1e7], [0]);
% plot(t,y)
% O1D = y(length(y),1);
% HOx = 1.5*10^-10*O1D*4.2*10^11+2*2.2*10^-10*O1D*1.68*10^12;

%% b
% kOHO3 = 1.7*10^-12*exp(-940/227);
% kHO2O = 3.0e-11*exp(200/227);
% kHO2O3 = 1.0e-14*exp(-490/227);
% O3 = 4.998*10^12;
% O = 5.586*10^8; %from problem 1
% OHHO2 = (kHO2O3*O3+kHO2O*O)/(kOHO3*O3);

%% c
% kOHHO2 = 4.8e-11*exp(250/227);
% HO2 = sqrt(7.9960e4/(kOHHO2*0.3419))
% HO = 0.3419*HO2

%% d
% [t,y] = ode15s(@HW3_O3_solve, [0 1e7], [0 0]);
% plot(t,y(:,1))
% O3SS = y(length(y),1)

%% e
% [t,y] = ode15s(@HW3_O3_solve, [0 1e7], [10^12 0 10^7 10^7]);
% plot(t,y(:,1))
% O3SS = y(length(y),1)
% figure (2)
% plot(t,y(:,2))
% figure (3)
% plot(t,y(:,3))
% figure (4)
% plot(t,y(:,4))

%% f
% kO3 = 0.0012;
% kO3O = 9.16*10^-16;
% kOHO3 = 1.7*10^-12*exp(-940/227);
% kHO2O3 = 1.0e-14*exp(-490/227);
% Closs = -kO3*y(130,1)-kO3O*y(130,1)*y(130,2);
% HOxloss = -kOHO3*y(130,3)*y(130,1)-kHO2O3*y(130,4)*y(130,1);

%% Problem 3 a
% jN2O = 5e-8;
% kN2Oa = 6.7e-11;
% kN2Ob = 4.9e-11;
% k4 = 3.2e-11;
% M = 4.2*10^17;
% jO3O1D = 2.7*10^-4;
% O3 = 4.998*10^12;
% 
% O1D = jO3O1D*O3/(k4*M);
% NO = 2*kN2Oa*O1D/(jN2O+(kN2Oa+kN2Ob)*O1D)*155;

%% d
% kNO2OH = 3.1915e-13;
% jHNO3 = 1.5e-5;
% kOHO3 = 1.7*10^-12*exp(-940/227);
% kHO2NO = 1.05e-11;
% kO3NO = 4.05e-15;
% kO3NO2 = 2.465e-18;
% jNO2 = 1.5e-2;
% HO2 = 3.7171e8;
% OH = 4.2779e7;
% O3 = 3.1557e12;
% HNO3 = kNO2OH/jHNO3*(kOHO3*OH*O3/(kHO2NO*HO2))*(kO3NO*O3+kHO2NO*HO2)/(jNO2+kO3NO2*O3)*OH

%% f
% [t,y] = ode15s(@HW3_NOx_solve, [0 1e7], [0 0 0 0 0 0 0]);
% plot(t,y(:,1))
% title('O3')
% figure (2)
% plot(t,y(:,2))
% title('O')
% figure (3)
% plot(t,y(:,3))
% title('OH')
% figure (4)
% plot(t,y(:,4))
% title('HO2')
% figure (5)
% plot(t,y(:,5))
% title('O(1D)')
% figure (6)
% plot(t,y(:,6))
% title('NO')
% figure (7)
% plot(t,y(:,7))
% title('NO2')

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

dO3 = -jO3*y(198,1) -kO3O*y(198,2)*y(198,1) -kOHO3*y(198,1)*y(198,3) -kHO2O3*y(198,1)*y(198,4) -kO3NO*y(198,6)*y(198,1) -kO3NO2*y(198,7)*y(198,1); %O3
tao=-y(198,1)/dO3/60