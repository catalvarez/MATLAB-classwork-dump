g = 9.8; %m/s2
rho = 1025; %kg/m3

P = waterDepth.*9.8.*1025./10000;

S = salinity;
T = temperature;
PR = 0;
PT = sw_ptmp(S,T,P,PR);
plot(PT,waterDepth)
set(gca,'YDir','reverse')
xlabel('Potential Temperature')
ylabel('Depth (m)')
title('Potential Temperature Profile')