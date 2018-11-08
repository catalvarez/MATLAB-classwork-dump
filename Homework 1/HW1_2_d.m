tspan = [0, 5e7];
c0 = [9e-5];

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 1);

[t, C] = ode15s(@(t,C) HW1_2_2(t,C), tspan, c0, options);

plot(t,C)
xlabel('Time [years]')
ylabel('Concentration [mol/kg]')