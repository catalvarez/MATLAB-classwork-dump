G = 6.67*10^(-11);
rhoc = 2700;
rhow = 1000;
rhom = 3300;
t = 30000;
h = t*(rhom-rhoc)/(rhom-rhow);

sigma1 = h*(rhoc-rhow);
sigma2 = (t-h)*(rhoc-rhom);
b1 = h/2;
b2 = h+(t-h)/2;

x = -100000:100000;

g_anom1 = 2*G*sigma1.*(pi/2.+atan(-x./b1));
g_anom2 = 2*G*sigma2.*(pi/2.+atan(-x./b2));
g_total = g_anom1 + g_anom2;

figure(1)
plot(x,g_anom1)
hold on
figure(2)
plot(x,g_anom2)
hold on
figure(3)
plot(x,g_total)
xlabel('x [m]')
ylabel('gravity anomaly')
hold on
