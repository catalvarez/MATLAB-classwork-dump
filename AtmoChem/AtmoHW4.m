M = 4.2*10^17;
VOC = 80E-9*M;
NO = 50E-9*M;
NO2 = 30E-9*M;

[t,y] = ode45(@HW4_O3_solve, [0 43200], [VOC 0 0 NO NO2 0]);
plot(t,y(:,1))
set(gca,'XTick',0:3600:43200,'XTickLabel',{'6:00','7:00','8:00','9:00','10:00'...
    '11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00'})
title('Initial Simulation')
ylabel('Concentration [molec/cm^3]')
xlabel('Time')
hold on
plot(t,y(:,4)+y(:,5),'g')
plot(t,y(:,6),'r')
legend('VOC','NOx','O3')