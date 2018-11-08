%   This script uses 'Definitions_ThreeBox_PO4.m' and 'Equations_ThreeBox_PO4.m' to
%   then solve for the time dependent solution.

%   Set the initial conditions

Definitions_ThreeBox_PO4_AtmCO2_2012;

%PO4
c0 = zeros(4,1);
c0(po4_1) = 2.5e-3;        % These [PO4] units are mol/m^3
c0(po4_2) = 2.5e-3;
c0(po4_3) = 2.5e-3;
c0(po4_4) = 0;

%options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 3.14e7);

c=ode15s(@(t,c) Equations_FourBox_PO4_2013(t,c),[0 16000*3.14e7],c0, options);
    
ConcMatrix_po4 = c.y(:,end)

%Alk
A0 = zeros(4,1);
A0(1) = 2375*10^-6;
A0(2) = 2375*10^-6;
A0(3) = 2375*10^-6;
A0(4) = 0;

%options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 3.14e7);

A=ode15s(@(t,A) Equations_FourBox_Alk_2013(t,A),[0 16000*3.14e7],A0, options);

ConcMatrix_Alk = A.y(:,end)

%DIC
D0 = zeros(4,1);
D0(1) = 2260*10^-6;
D0(2) = 2260*10^-6;
D0(3) = 2260*10^-6;
D0(4) = 0;

%options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 3.14e7);

D=ode15s(@(t,D) Equations_FourBox_DIC_2013(t,D),[0 16000*3.14e7],D0, options);

ConcMatrix_DIC = D.y(:,end)


global pco2
[co2,pco2,co3,ph,kh,o2,kspc,kspa] = dafunPE_jfa(ConcMatrix_DIC(2),ConcMatrix_Alk(2),TC(1),34,P(1))

%pCO2
p0 = zeros(4,1);
p0(1) = 0;
p0(2) = 0;
p0(3) = 0;
p0(4) = 280;

%options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', 3.14e7);

p=ode15s(@(t,p) Equations_FourBox_pCO2_2013(t,p),[0 16000*3.14e7],p0, options);

ConcMatrix_pCO2 = p.y(:,end)