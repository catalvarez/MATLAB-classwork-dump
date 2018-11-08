function [ dT] = Z( t,T )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dT = zeros(1,1);
sigma = 5.67.*10.^-8;
F0 = 1367;
T0 = (F0./(4.*sigma)).^(1/4);
A = 0.29;
tao = 0.63;
epsilon = 1/(1+tao);
Cp = 1000;
P = 10^5;
g = 10;
C = Cp.*P./g;

dT(1,1) = sigma/C.*(T0.^4.*(1-A)-epsilon.*T(1,1).^4).*86400;
% dT(1,1) = sigma/C.*(T0.^4.*(1-(0.47-0.25.*tanh((T(1,1)-268)./23)))-epsilon.*T(1,1).^4).*86400;


end

