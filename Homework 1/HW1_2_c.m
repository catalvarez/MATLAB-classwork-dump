tau = [1/5 1 5].*3.8e6;
phaseshift = [(9.499e6-9.859e6) (2.088e7-2.269e7) (1.609e8-1.709e8)];
A = [(1.903e-5-1.627e-5)/2 (9.523e-5-8.139e-5)/2 (0.0004762-0.0004071)/2];
plot(tau,phaseshift)
ylabel('Phase Shift')
xlabel('tau')