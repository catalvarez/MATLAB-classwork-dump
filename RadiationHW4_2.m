resol = [0.003 0.01 0.03 0.1 0.3];

v = 6239.9001:0.0001:6247.9;
P = logspace(6,5,100);
Sv0 = 1.53e-23; %cm
v0 = 6243.9; %cm-1
A = 7.2e-5/1000; %cm-1*(g*s-2*cm-1)-1
Ps = 1000*1000; %(g*s-2*cm-1)
m = 29*1.66e-24; %grams
g = 1000; %cm*s-2
X = 400/10^6; %mixing

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

H8 = zeros(5,100);
A8 = zeros(5,1);
for k = 1:5
%Interpolation
vi = 6239.901:resol(k):6247.9;
Ti = zeros(length(vi),1);
for i = 1:length(vi)
    Ti(i) = Tc(40000+(i-1)*(resol(k)/0.0001)+1);   
end

% Noise
SNR = 100;
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

K = zeros(100,2667);
for i = 1:length(P)
K(i,:) = -Sv0*X*2*A^2*Ps*P(i)./(2*pi*m*g*A*(A^2*P(i).^2+(v-v0).^2)).*((A^2*Ps^2+(v-v0).^2)./(v-v0).^2).^(-Sv0*X/(2*pi*m*g*A));
end
K = K.';

Sh = (K.'*Se'*K+Sa')';
H8(k,:) = 1/2*log(abs(Sh'*Sa));
A2 = Sh*K.'*Se'*K;
DOF(k) = trace(A2);
end