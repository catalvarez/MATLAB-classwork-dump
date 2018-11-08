NO20 = 2.45*10^12;
jNO2 = 9.971*10^-3;
t = 1:900;
NO2 = NO20.*exp(-jNO2.*t);
NO = NO20-NO20.*exp(-jNO2.*t);
O = NO20-NO20.*exp(-jNO2.*t);

plot(t,NO2)
hold on
plot(t,NO,'r')%y,m,c,r,g,b,w,k
% plot(t,O,'g')
title('Photolysis')
xlabel('Time (s)')
ylabel('Concentration (molecules cm-3')
