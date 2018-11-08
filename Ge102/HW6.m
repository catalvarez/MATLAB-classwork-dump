w0 = 446;
l = 1.8*10^6;
p = 3300;
g = -9.8;
B = 3*10^-15;
m = 10^21;
% B = l/(2*m);
% B = 1*10^-21;
t = 0:100:9320;
n = -2.2;

% w = w0*exp(B*p*g*l*t);
% w = ((1-n)*B*l*(p*g).^n.*t+w0.^(1-n)).^(1/(1-n));



% plot(t,w)
% hold on

age = 9320-[9320,9200,8800,8000,6000,0];
w2 = [446,437,385,331,261,180];

line = lsqcurvefit(HW6_trial(t),w0,age,w2);

scatter(age,w2)

