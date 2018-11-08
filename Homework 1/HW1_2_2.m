function [ dCdt] = HW1_2_2( t,C )
% 1-box model of the ocean


F0 = 3.3e10;
C0 = 90e-6;
A = 0.09.*F0;
%tau = 10^5;
tau = 3.8e6;
k = 1./tau;
M = 1.42e21;

dCdt = (F0 + A.*sin(2.*pi.*t./tau))./M - k.*C;

end

