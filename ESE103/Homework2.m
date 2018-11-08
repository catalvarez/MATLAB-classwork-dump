y0 = [80, 30, 30]
tspan = [1 150]

[t,y] = ode45(@(t,y) NPZmodel(t,y),tspan,y0)

semilogx(t,y)