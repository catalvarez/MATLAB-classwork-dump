%% Imports data
mat = dlmread('Europa_data.txt');
E_w = mat(:,1);
E_r  = mat(:,2);
E_sigma = mat(:,3);
mat = dlmread('snow1a.txt');
s_w = mat(:,1);
s_r  = mat(:,2);
s_sigma = mat(:,3);
mat = dlmread('H2O_ice_gds_77K.txt');
i_w = mat(:,1);
i_r  = mat(:,2);
i_sigma = mat(:,3);
figure (1)
errorbar(E_w,E_r,E_sigma,'.')
hold on
plot(s_w,s_r,'g','Linewidth',2)
plot(i_w,i_r,'r','Linewidth',2)
legend('Europa','Snow (278K)','Frost (77K)','Location',[0.6 0.75 0.2 0.1])
legend boxoff
title('Raw Spectral Comparison','FontSize',16)
xlabel('Wavelength (microns)','FontSize',14)
ylabel('Reflectance','FontSize',14)
%% Anaylsis
E_w_m1 = round(E_w.*1000)./1000;
L_m1 = length(E_w);
s_r_E = zeros(L_m1,1);
for i = 1:L_m1
    [~,I] = min(abs(s_w-E_w_m1(i)));
    s_r_E(i)=s_r(I);
end

a = 0:0.0001:0.5;
chi2 = zeros(length(a),1);
for j = 1:length(a)
chi2(j) = sum((s_r_E(26:101)-(E_r(26:101)-a(j))).^2./E_sigma(26:101).^2);
end
[x,I] = min(chi2);
t1 = a(I);
figure (2)
plot(E_w_m1(26:101),E_r(26:101),s_w(651:end),s_r(651:end),E_w_m1(26:101),s_r_E(26:101),'g+',E_w_m1(26:101),E_r(26:101),'b+')
title('Model 1 (Snow) No Correction','FontSize',16)
legend('Europa','Snow (278K)','Location',[0.6 0.75 0.2 0.1])
legend boxoff
xlabel('Wavelength (microns)','FontSize',14)
ylabel('Reflectance','FontSize',14)
figure (3)
plot(E_w_m1(26:101),E_r(26:101)-t1,s_w(651:end),s_r(651:end),E_w_m1(26:101),s_r_E(26:101),'g+',E_w_m1(26:101),E_r(26:101)-t1,'b+')
title('Model 1 (Snow) Best Fit','FontSize',16)
legend('Europa','Snow (278K)','Location',[0.6 0.75 0.2 0.1])
legend boxoff
xlabel('Wavelength (microns)','FontSize',14)
ylabel('Reflectance','FontSize',14)
figure (4)
plot(a,chi2,'Linewidth',3)
title('Chi-squared Minimization, Model 1','FontSize',16)
xlabel('Subtraction Factor, b','FontSize',14)
ylabel('Chi-squared Value','FontSize',14)

%% Model 2
E_w_m2 = round(E_w.*1000)./1000;
L_m2 = length(E_w);
i_r_E = zeros(L_m2,1);
for i = 1:L_m2
    [~,I] = min(abs(i_w-E_w_m2(i)));
    i_r_E(i)=i_r(I);
end

b = 0:0.0001:0.4;
chi2_2 = zeros(length(b),1);
for i = 1:length(b)
chi2_2(i) = sum((i_r_E(26:101)-(E_r(26:101)-b(i))).^2./E_sigma(26:101).^2);
end
[x2,I] = min(chi2_2);
t2 = b(I);
figure (5)
plot(E_w_m2(26:101),E_r(26:101),i_w(18:225),i_r(18:225),'r',E_w_m2(26:101),i_r_E(26:101),'r+',E_w_m2(26:101),E_r(26:101),'b+')
title('Model 2 (Frost) No Correction','FontSize',16)
legend('Europa','Frost (77K)','Location',[0.6 0.75 0.2 0.1])
legend boxoff
xlabel('Wavelength (microns)','FontSize',14)
ylabel('Reflectance','FontSize',14)
figure (6)
plot(E_w_m2(26:101),E_r(26:101)-t2,i_w(18:225),i_r(18:225),'r',E_w_m2(26:101),i_r_E(26:101),'r+',E_w_m2(26:101),E_r(26:101)-t2,'b+')
title('Model 2 (Frost) Best Fit','FontSize',16)
legend('Europa','Frost (77K)','Location',[0.6 0.75 0.2 0.1])
legend boxoff
xlabel('Wavelength (microns)','FontSize',14)
ylabel('Reflectance','FontSize',14)
figure (7)
plot(b,chi2_2,'Linewidth',3)
title('Chi-squared Minimization, Model 2','FontSize',16)
xlabel('Subtraction Factor, b','FontSize',14)
ylabel('Chi-squared Value','FontSize',14)
