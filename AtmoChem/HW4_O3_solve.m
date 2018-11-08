function dy = HW4_O3_solve(t,y)
% Simple 0D model of tropospheric ozone throughout the day
% y = [RH RO2 HO2 NO NO2 O3]

OH = 5E6;
kRHOH = 26.3E-12;
kRO2NO = 7.7E-12;
kHO2NO = 8.1E-12;
kOHNO2 = 1.1E-11;
kHO2HO2 = 2.9E-12;
kRO2HO2 = 5.2E-12;
kO3NO = 1.9E-14;
% photolysis of NO2 throughout the day
i = floor(t./3600);
if i == 0
    jNO2 = 2.215E-3;
end
if i == 1
    jNO2 = 4.645E-3;
end
if i == 2
    jNO2 = 6.646E-3;
end
if i == 3
    jNO2 = 8.108E-3;
end
if i == 4
    jNO2 = 9.088E-3;
end
if i == 5
    jNO2 = 9.634E-3;
end
if i == 6
    jNO2 = 9.766E-3;
end
if i == 7
    jNO2 = 9.492E-3;
end
if i == 8
    jNO2 = 8.799E-3;
end
if i == 9
    jNO2 = 7.659E-3;
end
if i == 10
    jNO2 = 6.014E-3;
end
if i == 11
    jNO2 = 3.831E-3;
end
if i == 12
    jNO2 = 1.421E-3;
end

r1 = kRHOH*y(1)*OH;
r2 = kRO2NO*y(2)*y(4);
r3 = kHO2NO*y(3)*y(4);
r4 = kOHNO2*OH*y(5);
r5 = kHO2HO2*y(3)*y(3);
r6 = kRO2HO2*y(2)*y(3);
r7 = jNO2*y(5);
r8 = kO3NO*y(6)*y(4);
dy(1,1)= -r1;%RH
dy(2,1)= r1 -r2 -r6;%RO2
dy(3,1)= r2 -r3 -2*r5 -r6;%HO2
dy(4,1)= -r2 -r3 +r7 -r8;%NO
dy(5,1)= r2 +r3 -r4 -r7 +r8;%NO2
dy(6,1)= r7 -r8;%O3

end