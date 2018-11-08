%   This script solves simple parent-daughter radioactive decay series
%   using ODE solvers and plots them with the analytical solutions.  


clear all 
close all
clc

% Set the initial conditions

decaydefs; %call the definitions file 

c0 = zeros(2,1);

c0(N) = 5e3; %Start with lots of N and no M
c0(M) = 0;

tinit = 0;
tfinal = 100;
tstep = 1;

% Set up the ODE.  Use ODE 15s.  

options = odeset('RelTol', 1e-5, 'AbsTol', 1e-10, 'InitialStep', tstep);

c=ode15s(@(t,c) parentdaughterdecay(t,c),[tinit tfinal],c0, options);

t = c.x;
N_t = c.y(N,:);
M_t = c.y(M,:);

% Write out analytical solutions for decay of N and M.  

Na_t = c0(N).*exp(-t./tao_n);
Ma_t = c0(N).*tao_m./(tao_n-tao_m).*(exp(-t./tao_n)-exp(-t./tao_m));

%Find the difference at each time step

Ndiff = Na_t - N_t;
Mdiff = Ma_t - M_t;

%plot the data. 

figure()
subplot(2,1,1)
plot(t,N_t,'ro')
hold on
plot(t,M_t,'o')
plot(t,Na_t,'g')
plot(t,Ma_t,'k')
hold off
xlabel('time')
ylabel('Amount')
legend('N','M','analytical N','analytical M')
title('Time evolution of N and M')

subplot(2,1,2)
plot(t,Ndiff,'r')
hold on
plot(t,Mdiff)
xlabel('time')
ylabel('difference')
title('Difference between numerical and analytical solutions')
