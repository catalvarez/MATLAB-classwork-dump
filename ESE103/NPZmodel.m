function dy = NPZmodel(t,y)
% NPZmodel: Ocean nitrification model
% t: time vector, in units of day of year
% y: initial conditions for N, P and Z (as a vector)
% dy returns a matrix of Y values of size (length(timevec),3), 
% where Y(:,1) is N, Y(:,2) is P, and Y(:,3) is Z. 

% Assign arguments %
N = y(1); % Nitrogen
P = y(2); % Phytoplankton
Z = y(3); % Zooplankton

Vmax =1; % default = 1
k_N = 20;
k_ZP = 0.01; % default = 0.01
k_Z = 0.3;
lambda_Z = 0.5;

dy = zeros(3,1); 
dy(1) = (lambda_Z*Z + k_ZP*k_Z*P*Z) - Vmax*N/(N+k_N)*P;
dy(2) = Vmax*N/(N+k_N)*P - k_ZP*P*Z;
dy(3) = k_ZP*(1-k_Z)*P*Z - lambda_Z*Z;

%{
% Calculating light limitations 
%(you need to put this in the correct place in the NPZ equations above)
doy = floor(mod(t,364)); % Converts t to day of year
lat = 30; % latitude in degrees (default is 30)
S = daily_insolation(doy,1367,lat); % Determines local insolation 
Sm = 500;
L = exp(1/Sm - 1/S); %determines the coefficient for modifying growth under light limitation
%}