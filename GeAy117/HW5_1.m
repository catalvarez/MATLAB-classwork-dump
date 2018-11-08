[a,~,~,~,~,H] = fminunc(@StatsHW5,[1 4 4]);

covar = inv(H);

rho0_sigma = sqrt(covar(1,1));
K_prime_sigma = sqrt(covar(2,2));
K0_sigma = sqrt(covar(3,3));
mat = dlmread('ps5_data2.txt');

figure (1)
line = a(1).*(1+a(2).*mat(:,1)./a(3)).^(1./a(2));
errorbar(mat(:,1),mat(:,2),mat(:,3),'.')
hold on
plot(mat(:,1),line,'g')
title('Optimal Estimation')
xlabel('P (GPa)')
ylabel('density (g/cm^3)')

correl = [1 covar(1,2)/(rho0_sigma*K_prime_sigma) covar(1,3)/(rho0_sigma*K0_sigma);...
    covar(2,1)/(rho0_sigma*K_prime_sigma) 1 covar(2,3)/(K0_sigma*K_prime_sigma);...
    covar(3,1)/(rho0_sigma*K0_sigma) covar(3,2)/(K0_sigma*K_prime_sigma) 1];

bound95 = chi2inv(.95,2);
bound68 = chi2inv(.68,2);

mat = dlmread('ps5_data2.txt');
grid_rho = 3.32:0.001:3.36;
grid_kp = 2:0.1:7;
grid_k0 = 115:150;
dchi2_rho_kp = zeros(length(grid_rho),length(grid_kp));
dchi2_rho_k0 = zeros(length(grid_rho),length(grid_k0));
dchi2_kp_k0 = zeros(length(grid_kp),length(grid_k0));
for j = 1:length(grid_kp)
    for k = 1:length(grid_k0)
        for i = 1:length(grid_rho)
            b1 = [grid_rho(i)-a(1); grid_kp(j)-a(2)];
            covarmarg1 = [covar(1,1) covar(1,2); covar(2,1) covar(2,2)];
            dchi2_rho_kp(i,j) = abs(b1.'*(covarmarg1^(-1))*b1);
            
            b2 = [grid_rho(i)-a(1); grid_k0(k)-a(3)];
            covarmarg2 = [covar(1,1) covar(1,3); covar(3,1) covar(3,3)];
            dchi2_rho_k0(i,k) = abs(b2.'*(covarmarg2^(-1))*b2);
        
            b3 = [grid_kp(j)-a(2); grid_k0(k)-a(3)];
            covarmarg3 = [covar(2,2) covar(2,3); covar(3,2) covar(3,3)];
            dchi2_kp_k0(j,k) = abs(b3.'*(covarmarg3^(-1))*b3);
    % input a = rho0,K_prime,K0
        end
    end
end
figure (2)
subplot(2,2,1)
contour(grid_kp,grid_rho.',dchi2_rho_kp,[bound68, bound95])
title('Density vs dK/dP')
xlabel('K prime')
ylabel('rho')
hold on
plot(a(2),a(1),'x')
subplot(2,2,2)
contour(grid_k0,grid_rho.',dchi2_rho_k0,[bound68, bound95])
title('Density vs Bulk Modulus')
xlabel('K_0')
ylabel('rho')
hold on
plot(a(3),a(1),'x')
subplot(2,2,3)
contour(grid_k0,grid_kp.',dchi2_kp_k0,[bound68, bound95])
title('dK/dP vs Bulk Modulus')
xlabel('K_0')
ylabel('K prime')
hold on
plot(a(3),a(2),'x')

