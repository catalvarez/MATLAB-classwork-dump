function [H, DOF] = RadiationHW4_3(resol,SNR)
v = 6239.9001:0.0001:6247.9;
P = logspace(6,5,100);
Sv0 = 1.53e-23; %cm
v0 = 6243.9; %cm-1
A = 7.2e-5/1000; %cm-1*(g*s-2*cm-1)-1
Ps = 1000*1000; %(g*s-2*cm-1)
m = 29*1.66e-24; %grams
g = 1000; %cm*s-2
X = 400/10^6; %mixing
%resol = 0.01;

%Transmittance
T = ((A^2*Ps^2+(v-v0).^2)./(v-v0).^2).^(-Sv0*X/(2*pi*m*g*A));

% Instrumental Profile
FWHM = 0.01; %cm^-1
mean = 6243.9; %cm^-1
sig = (FWHM/2)^2/(2*log(2));
Gaus = 1/(sig*sqrt(2*pi))*exp(-(v-mean).^2/(2*sig^2));
Gaus = Gaus/sum(Gaus);
% Convolution
Tc = conv(Gaus,T);
%Interpolation
vi = 6239.9+resol:resol:6247.9;
Ti = zeros(length(vi),1);
for i = 1:length(vi)
    Ti(i) = Tc(40000+(i-1)*(resol/0.0001)+1);   
end

% Noise
%SNR = 100;
error = Ti/SNR;
Se = zeros(length(Ti));
sig = error.^2;
for i = 1:length(Ti)
    Se(i,i) = sig(i);
end

% Prior Uncertainty
Sa = zeros(100);
sig = (0.03*X)^2;
for i = 1:100
    Sa(i,i) = sig;
end

K = zeros(100,80000);
for i = 1:length(P)
K(i,:) = -Sv0*X*2*A^2*Ps*P(i)./(2*pi*m*g*A*(A^2*P(i).^2+(v-v0).^2)).*((A^2*Ps^2+(v-v0).^2)./(v-v0).^2).^(-Sv0*X/(2*pi*m*g*A));
end
K = K.';

%Conv and Interp
Ki = zeros(length(vi),length(P));
Ti = zeros(length(vi),1);
for i = 1:length(P)
Tc = conv(Gaus,K(:,i));
for j = 1:length(vi)
    Ti(j) = Tc(40000+(j-1)*(resol/0.0001)+1);   
end
Ki(:,i) = Ti;
end

Sh = (Ki.'*Se'*Ki+Sa')';
H = 1/2*log(det(Sh'*Sa));
A2 = Sh*Ki.'*Se'*Ki;
DOF = trace(A2);