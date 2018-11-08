[a,~,~,~,~,H] = fminunc(@StatsHW4,[1 0 124000]);

covar = inv(H);

msigma = sqrt(covar(1,1));
csigma = sqrt(covar(2,2));
Psigma = sqrt(covar(3,3));
mat = dlmread('ps4_data.txt');

tline = 0:3720:372000;
line = a(1).*sin(2.*pi./a(3).*tline)+a(2);
errorbar(mat(:,1),mat(:,2),mat(:,3),'o')
hold on
plot(tline,line,'g')
title('Method Two: Numerical Best Fit')
xlabel('time (years)')
ylabel('obliquity (degrees)')
