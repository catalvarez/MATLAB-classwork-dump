function dy = ps3prob2e(t,y)

N2O= 65100000000;
j1 = 5*10^-8;
j2 = 1.5*10^-5;
j3 = .00027;
j4 = 1.5*10^-2;
k1 = 6*10^-11;
k2 = 6*10^-34;
k3 = .0012;
k4 = 9.16*10^-16;
k5 = 1.7*10^-12*exp(-940/227);
k6 = 1.0*10^-14*exp(-490/227);
k7 = 3.0*10^-11*exp(200/227);
k8 = 4.8*10^-11*exp(250/227);
k9 = 3.0*10^-12*exp(-1500/227);
k10 = 1.2*10^-13*exp(-2450/227);
k11 =5.6*10^-12*exp(180/227);
k12 = 2.2*10^-10;
k13 = 3.5*10^-12*exp(210/227);
k14 = 6.7*10^-11;
k15 = 4.9*106-11;
k16 = 1.5*10^-10;
k17 = 3.2*10^-11*exp(70/227);
k18 = 1.8*10^-11*exp(110/227);
k19 = 3.1915*10^-13;
M = 4.2*10^17;
O2 = .21*M;
N2 = M*.79;
H2O = 1.68*10^12;
HNO3 = 2.5113e+09;
CH4 = 4.2*10^11;
dy = zeros(7,1);

dy(1) = k2*O2*y(2)*M - k3*y(1) - k4*y(2)*y(1) - k5*y(3)*y(1) - k6*y(4)*y(1) - k10*y(7)*y(1) - k9*y(5)*y(1); %do3

dy(2) = 2*k1*O2 + k3*y(1) + j4*y(7) - k4*y(2)*y(1) - k2*O2*y(2)*M - k7*y(2)*y(4) - k11*y(7)*y(2); %do

dy(3) = k7*y(4)*y(2) - k5*y(3)*y(1) - k19*y(7)*y(3) + k6*y(4)*y(1) - k8*y(3)*y(4) + 2*k12*y(6)*H2O + j2*HNO3 + k13*y(4)*y(5) +k16*y(6)*CH4; %dOH

dy(4) = k5*y(3)*y(1) - k6*y(4)*y(1) - k8*y(3)*y(4) - k13*y(4)*y(5); %dHO2

dy(5) = 2*k14*y(6)*N2O - k9*y(5)*y(1) + j4*y(7) + k11*y(7)*y(2)- k13*y(4)*y(5);%dNO

dy(6) = j3*y(1) - k16*y(6)*CH4 - k12*y(6)*H2O - k17*y(6)*O2 - k18*y(6)*N2 + j1*N2O - k15*y(6)*N2O - k14*y(6)*N2O; %dO1

dy(7) = j2*HNO3 - j4*y(7)- k19*y(7)*y(3) - k11*y(7)*y(2) - k10*y(7)*y(1) +k9*y(5)*y(1) + k13*y(4)*y(5);%dNO2
end