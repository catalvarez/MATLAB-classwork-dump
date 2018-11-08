% Homework one: climate feedback

[t,T] = ode45('dT',[0 365],[288]);

plot(t,T)
title('Feedback Off')
axis([0 400 288 290])
xlabel('Time [days]')
ylabel('T [K]')
