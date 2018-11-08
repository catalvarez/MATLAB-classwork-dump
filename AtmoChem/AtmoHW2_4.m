tspan = [0 3000];
% NO2 = [2.45*10^12 0 0 0];
NO2 = [1.5485*10^12 0 0.9015*10^12 0.9015*10^12];

%options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 1);

[t, y] = ode15s(@AtmoHW2_4_solve, tspan, NO2); %, options);

plot(t,y)
title('Removed from Sunlight')
xlabel('Time [s]')
ylabel('Concentration [molecules cm-3]')