%   This script uses 'Definitions_ThreeBox_PO4.m' and 'Equations_ThreeBox_PO4.m' to
%   then solve for the time dependent solution.

%   Set the initial conditions

Definitions_ThreeBox_PO4_2013;

c0 = zeros(3,1);

c0(po4_1) = 2.5e-3;        % These [PO4] units are mol/m^3
c0(po4_2) = 2.5e-3;
c0(po4_3) = 2.5e-3;

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 3.14e7);

c=ode15s(@(t,c) Equations_ThreeBox_PO4_2013(t,c),[0 16000*3.14e7],c0, options);
    
ConcMatrix = c.y(:,end)

%Figures_ThreeBox_PO4;