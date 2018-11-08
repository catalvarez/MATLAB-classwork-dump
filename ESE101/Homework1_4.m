% Homework one: climate feedback
sigma = 5.67.*10.^-8;
F0 = 1367;
%T0 = (F0./(4.*sigma)).^(1/4);
A = 0.29;
%tao = 0.63;
%epsilon = 1/(1+tao);
Cp = 1000;
P = 10^5;
g = 10;
C = Cp.*P./g;
T0 = 260:1:360;
    F1 = T0.^4.*(1-A);
T = 260:1:360;
    F2 = T.^4;
    
Z = ones(101,101);
for i = 1 : 101
    for j = 1 : 101
        Z(i,j) = F1(j)-F2(i)./(1+(0.56+0.07.*exp(5413.*(1/288-1./T(i)))));
    end
end
%Z = sigma/C.*(T0.^4.*(1-A)-epsilon.*T.^4);
Y = sigma/C.*Z.*60.*60.*24;
[M,h] = contour(T0,T,Y)
clabel(M,h)
xlabel('T_0 [K]')
ylabel('T [K]')