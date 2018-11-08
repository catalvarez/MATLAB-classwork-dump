% Ge102 HW2 Winter 2014
%
% Find the radial density structure of earth by integrating
% Adams-Williamson equation

[r, z, a, b] = textread('PREM.txt', '%f %f %f %f', ...
    'headerlines', 5, 'commentstyle', 'shell');

% Remaining mass, density and gravity at different depths
M_rem = zeros(1, length(r));
M_shell = zeros(1, length(r));
d = zeros(1, length(r));
g = zeros(1, length(r));
I = zeros(1, length(r));
r = r*1000;
z = z*1000;

Me          = 5.974e24; % Mass of earth, in kg
M_rem(1)    = Me;    % Remaining mass of earth
d(1)        = 3300;  % surface density, in kg/m^3
g(1)        = 9.81;  % surface gravity, in m/s^2
inertia     = 0.0;

% Loop over depth segment from surface down to core
for i = 2:length(r)
    % TODO: Calculate seismic parameter, change in radius, and hence
    %       the shell mass, remaining mass, new gravity and new density
    % sp = seismic parameter
    sp = a(i)^2-4/3*b(i)^2;
    dr = r(i-1)-r(i);
    M_shell(i) = 4/3*pi*(r(i-1)^3-r(i)^3)*d(i-1);
    M_rem(i) = M_rem(i-1)-M_shell(i);
%     g(i) = 6.67*10^-11*M_rem(i)/(r(i))^2;
    g(i) = 6.67*10^-11*M_rem(i)/(6370000)^2;
%  g(i) = g(i-1)-6.67*10^-11*(M_shell)/(r(i-1)-r(i))^2;
    d(i) = d(i-1)/(1+g(i)/sp*dr);
    
    if r(i) == 5970821
        d(i) = d(i)*1.05;
    end
    if r(i) == 5701911
        d(i) = d(i)*1.09;
    end
    if r(i) == 3480536
        d(i) = 9900;
    end
    
    
    % TODO: Account for density jumps at ICB, CMB, 670 and 410
end

% Make plots
plot(r,a)
xlabel('Radius (m)')
ylabel('a (m/s)')
figure
plot(r,b)
xlabel('Radius (m)')
ylabel('b (m/s)')
figure
plot(r,g)
xlabel('Radius (m)')
ylabel('g (m/s^2)')
figure
plot(r,d)
xlabel('Radius (m)')
ylabel('density (kg/m^3)')
