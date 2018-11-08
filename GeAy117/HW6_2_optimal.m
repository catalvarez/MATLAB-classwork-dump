[a,~,~,~,~,H] = fminunc(@HW6_2,[70 .85 -26]);

covar = inv(H);

K_sigma = sqrt(covar(1,1));
phi_sigma = sqrt(covar(2,2));
gamma_sigma = sqrt(covar(3,3));
mat = dlmread('ps6_rv_data.txt');

figure (1)
line = a(1).*cos(2*pi/6.5.*mat(:,1)+a(2))+a(3);
errorbar(mat(:,1),mat(:,2),mat(:,3),'.')
hold on
plot(mat(:,1),line,'g')
title('Optimal Estimation')
xlabel('P (GPa)')
ylabel('density (g/cm^3)')