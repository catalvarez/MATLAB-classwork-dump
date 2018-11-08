% Homework one: climate feedback
sigma = 5.67.*10.^-8;
F0 = 1367;
%T0 = (F0./(4.*sigma)).^(1/4);
%A = 0.30;
tao = 0.63;
epsilon = 1/(1+tao);
Cp = 1000;
P = 10^5;
g = 10;
C = Cp.*P./g;
T0 = 220:1:320;
    F1 = T0.^4;
T = 220:1:320;
    F2 = epsilon.*T.^4;    
    
Z = ones(101,101);
for i = 1 : 101
    for j = 1 : 101
        Z(i,j) = F1(j).*(1-(0.47-0.25.*tanh((T(i)-268)./23)))-F2(i);
    end
end
%Z = tao/C.*(T0.^4.*(1-A)-epsilon.*T.^4);
Y = sigma/C.*Z.*60.*60.*24;
[M,h] = contour(T0,T,Y)
clabel(M,h)
xlabel('T_0 [K]')
ylabel('T [K]')