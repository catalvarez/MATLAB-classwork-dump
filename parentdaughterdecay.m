function [ dcdt ] = parentdaughterdecay( t,c )
%A function to calculate the decay equation using an ode.

decaydefs; %call the definitions file

%Set up dcdt matrix for both n and m

dcdt(1) = - c(1)./tao_n;

dcdt(2) = c(1)./tao_n - c(2)./tao_m;

dcdt = dcdt';

end

