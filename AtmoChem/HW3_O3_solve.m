function dy = HW3_O3_solve(t,y)

kO2 = 6.0*10^-11;
kO2OM = 6.0*10^-34;
kO3 = 0.0012;
kO3O = 9.16*10^-16;
M = 4.2*10^17;
O2 = 0.21*M;
kOHO3 = 1.7*10^-12*exp(-940/227);
kHO2O3 = 1.0e-14*exp(-490/227);
% y(3) = 4.0245e7; % part 2d only
% y(4) = 1.3760e7; % part 2d only
kHO2O = 3.0e-11*exp(200/227);
kHO2OH = 4.8e-11*exp(250/227);

dy(1,1) = kO2OM*O2*y(2)*M-kO3*y(1)-kO3O*y(2)*y(1)-kOHO3*y(1)*y(3)-kHO2O3*y(1)*y(4); %O3
dy(2,1) = 2*kO2*O2+kO3*y(1)-kO3O*y(2)*y(1)-kO2OM*O2*y(2)*M-kHO2O*y(4)*y(2); %O
dy(3,1) = -kOHO3*y(3)*y(1)+kHO2O*y(4)*y(2)+kHO2O3*y(4)*y(1)-kHO2OH*y(3)*y(4);  %OH
dy(4,1) = -kHO2OH*y(3)*y(4)+kOHO3*y(3)*y(1)-kHO2O3*y(4)*y(1);   %HO2

end
