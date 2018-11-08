sig = 0.1:1000;
h = 6.626*10^-34;
c = 2.998*10^8;
k = 1.38*10^-23;
T = 280;

%for i = drange(0.1:100)
test = 1./(exp(h.*c./(sig.*k.*T))-1);

B = 2.*h.*c.^2./(sig.^5.*(exp(h.*c./(sig.*k.*T))-1));
B1 = B/max(B);

plot(sig,test)
figure (2)
semilogy(sig,B1)
