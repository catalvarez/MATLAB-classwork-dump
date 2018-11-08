g = 9.8; %m/s2
rho = 1025; %kg/m3

P = waterDepth.*9.8.*1025./10000;

S = salinity;
T = temperature;
PR = 0;
D = sw_pden(S,T,P,PR);
PT = sw_ptmp(S,T,P,PR);
plot(S,PT)
%set(gca,'YDir','reverse')
xlabel('Salinity')
ylabel('Potential Temperature')
title('Water Masses')