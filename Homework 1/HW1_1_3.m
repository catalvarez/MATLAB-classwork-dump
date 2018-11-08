g = 9.8; %m/s2
rho = 1025; %kg/m3

P = waterDepth.*9.8.*1025./10000;

S = salinity;
T = temperature;
PR = 0;
D = sw_pden(S,T,P,PR);
plot(D-2000,waterDepth)
set(gca,'YDir','reverse')
xlabel('Sigma-2')
ylabel('Depth (m)')
title('Sigma-2 Profile')