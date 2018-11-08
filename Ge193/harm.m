function [phase,ampl]=harm(obs,w,t);
% finds the amplitude and the phase of a sin curve of a given angular
% freq w that  gives best least square fit to obs collected at t

% matrix of basis functions
s=[sin(w*t(:)) cos(w*t(:))];

% solve for coefficients
coeff=s\obs(:);
ampl=norm(coeff);
phase=atan2(coeff(2),coeff(1))*180/pi;