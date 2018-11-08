% % IMPORTANT: remember to modify log_posterior for correct y=f(x)
% param = [71.8703 6.5 0.3473 23.4585]; %initial guesses
% stepsize = [0.2938 0.0001 0.0032 0.1869]; %step-sizes
% s = 100000; %number of steps per chain
% %Imports data
% mat = dlmread('ps6_rv_data.txt');
% X = mat(:,1)-mat(1,1);
% Y  = mat(:,2);
% sigma = mat(:,3);
% chain = zeros(s,length(param));
% accepts = 0;
% rejects = 0;
% 
% for j = 1:s
%     priorL = log_posterior(param,X,Y,sigma);
%     index = perturb_param(param);
%     trial = new_param(param,stepsize,index);
%     postL = log_posterior(trial,X,Y,sigma);
%     eval = evaluate_step(priorL,postL);
%     if eval == 1
%         chain(j,:) = trial;
%         param = trial;
%         accepts = accepts+1;
%     else
%         chain(j,:) = param;
%         rejects = rejects+1;
%     end
%     
% end
% 
% chainc_1 = chain(3000:end,1);
% chainc_2 = chain(2500:end,2);
% chainc_3 = chain(2500:end,3);
% chainc_4 = chain(3000:end,4);

% t = [mean(chainc_1) mean(chainc_2) mean(chainc_3) mean(chainc_4)];
% figure (1)
% hist(chainc_1,100)
% title('Velocity Semi-Amplitude Histogram')
% xlabel('Amplitude Value')
% ylabel('Counts')
% figure (2)
% plot(1:s,chain(:,1))
% title('Velocity Semi-Amplitude Trace')
% xlabel('Step')
% ylabel('K Value')
% K_conf = sqrt(sum((chainc_1-t(1)).^2)./s);
% figure (3)
% hist(chainc_2,100)
% title('Period Histogram')
% xlabel('P Value')
% ylabel('Counts')
% figure (4)
% P_conf = sqrt(sum((chainc_2-t(2)).^2)./s);
% plot(1:s,chain(:,2))
% title('Period Trace')
% xlabel('Step')
% ylabel('P Value')
% figure (5)
% hist(chainc_3,100)
% title('Phase Shift Histogram')
% xlabel('Phi Value')
% ylabel('Counts')
% figure (6)
% phi_conf = sqrt(sum((chainc_3-t(3)).^2)./s);
% plot(1:s,chain(:,3))
% title('Phase Shift Trace')
% xlabel('Step')
% ylabel('Phi Value')
% figure (7)
% hist(chainc_4,100)
% title('Vertical Offset Histogram')
% xlabel('Gamma Value')
% ylabel('Counts')
% figure (8)
% gamma_conf = sqrt(sum((chainc_4-t(4)).^2)./s);
% plot(1:s,chain(:,4))
% title('Vertical Offset Trace')
% xlabel('Step')
% ylabel('Gamma Value')

% acceptance = accepts/(accepts+rejects);
% binsize1 = 0.002;
% values1 = 74:binsize1:76.5;
% bin1 = histc(chainc_1,values1);
% binsize2 = 0.00002;
% values2 = 6.4947:binsize2:6.4952;
% bin2 = histc(chainc_2,values2);
% binsize3 = 0.002;
% values3 = .84:binsize3:0.92;
% bin3 = histc(chainc_3,values3);
% binsize4 = 0.02;
% values4 = -26.6:binsize4:-25.2;
% bin4 = histc(chainc_4,values4);
% KP_post2D = zeros(length(bin1),length(bin2));
% Kphi_post2D = zeros(length(bin1),length(bin3));
% Kg_post2D = zeros(length(bin1),length(bin4));
% Pphi_post2D = zeros(length(bin2),length(bin3));
% Pg_post2D = zeros(length(bin2),length(bin4));
% phig_post2D = zeros(length(bin3),length(bin4));
% for i = 1:length(bin1)
%     for j = 1:length(bin2)
%         for k = 1:length(bin3)
%             for l = 1:length(bin4)
% KP_post2D(i,j) = bin1(i).*bin2(j);
% Kphi_post2D(i,k) = bin1(i).*bin3(k);
% Kg_post2D(i,l) = bin1(i).*bin4(l);
% Pphi_post2D(j,k) = bin2(j).*bin3(k);
% Pg_post2D(j,l) = bin2(j).*bin4(l);
% phig_post2D(k,l) = bin3(k).*bin4(l);
%             end
%         end
%     end
% end

phig_post2D = phig_post2D./sum(sum(phig_post2D*binsize3*binsize4));
contour(values4,values3,phig_post2D)
title('2D Posterior')
xlabel('phi')
ylabel('gamma')

% thing = [chainc_1 chainc_2(501:end)];
%     hist3(thing)

% fold = mod(X,t(2));
% x = 0:0.1:t(2);
% line = t(1).*cos(2*pi./t(2).*x+t(3))+t(4);
% errorbar(fold,Y,sigma,'.')
% xlim([0 2*pi])
% hold on
% plot(x,line,'g')
% title('MCMC Best Fit')
% xlabel('Period')
% ylabel('Radial Velocity')

% G = 6.67*10^-11;
% M = 1.70;
% Msini = (chainc_2(501:end)*24*3600./(2*pi*G)).^(1/3).*chainc_1.*M^(2/3);
% t(5) = mean(Msini);
% plot(1:length(Msini),Msini)
% Msini_conf = sqrt(sum((Msini-t(5)).^2)./s);