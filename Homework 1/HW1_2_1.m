%F0 = 3.3.*10.^10;
%C0 = 90.*10.^(-6);
%A = 10.^10;
%tau = 3.8.*10.^9;
%k = F0/C0;
tspan = [0, 5e8];
c0 = [9e-5];

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 1);

[t, C] = ode15s(@(t,C) HW1_2_2(t,C), tspan, c0, options);

plot(t,C)
title('Tao 5x Residence Time')
xlabel('Time [years]')
ylabel('Concentration [mol/kg]')