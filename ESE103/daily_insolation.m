function Sd = daily_insolation(doy, S0, phi)
% DAILY_INSOLATION: Daily-mean insolation at top of atmosphere.
% 
%   DAILY_INSOLATION(DAY, S0, PHI) returns the daily-mean insolation
%   at day of year DAY [1,...,364] for solar constant S0 and latitude
%   PHI [angle in radians].
  
delta = declination_angle(doy);
phi   = phi*pi/180;  
h0    = real(acos(-tan(phi).*tan(delta)));
Sd    = S0/pi * (h0.*sin(phi).*sin(delta) + cos(phi).*cos(delta).*sin(h0));

function delta = declination_angle(doy)

% declination angle [Hartmann, Appendix A]
  theta_d = 2*pi*doy/365;
  delta   =   .006918 ...
	    - .399912 * cos(theta_d)   + .070257 * sin(theta_d) ...
	    - .006758 * cos(2*theta_d) + .000907 * sin(2*theta_d) ...
	    - .002697 * cos(3*theta_d) + .001480 * sin(3*theta_d);
return
