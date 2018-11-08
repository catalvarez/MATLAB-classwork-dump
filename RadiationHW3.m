Sv0 = 1.53e-23; %cm
v0 = 6243.9; %cm-1
A = 7.2e-5/1000; %cm-1*(g*s-2*cm-1)-1
Ps = 1000*1000; %(g*s-2*cm-1)
m = 29*1.66e-24; %grams
g = 1000; %cm*s-2
X = 400/10^6; %mixing
v = 6240:0.1:6247;

T = ((A^2*Ps^2+(v-v0).^2)./(v-v0).^2).^(-Sv0*X/(2*pi*m*g*A));
figure (1)
plot(v,T)
xlabel('frequency')
ylabel('Transmittance')

figure (2)
hold on
P = 0:10000:1000000;
for v = 6243.6:.01:6244
J = -Sv0*X*2*A^2*Ps*P./(2*pi*m*g*A*(A^2*P.^2+(v-v0)^2)).*((A^2*Ps^2+(v-v0).^2)./(v-v0).^2).^(-Sv0*X/(2*pi*m*g*A));
plot(P,J)
xlabel('pressure (dyn cm^-2)')
ylabel('Jacobian')
end