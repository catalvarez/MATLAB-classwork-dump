param = [1.0 4.0]; %initial guesses
stepsize = [0.5 1]; %step-sizes
s = 10000; %number of steps per chain
%Imports data
mat = dlmread('ps6_line_data.txt');
X = mat(:,1);
Y  = mat(:,2);
sigma = mat(:,3);
chain = zeros(s,length(param));
accepts = 0;
rejects = 0;

for j = 1:s
    priorL = log_posterior(param,X,Y,sigma);
    index = perturb_param(param);
    trial = new_param(param,stepsize,index);
    postL = log_posterior(trial,X,Y,sigma);
    eval = evaluate_step(priorL,postL);
    if eval == 1
        chain(j,:) = trial;
        param = trial;
        accepts = accepts+1;
    else
        chain(j,:) = param;
        rejects = rejects+1;
    end
    
end
% 
% t = mean(chain,1);
% figure (1)
% hist(chain(:,1),100)
% title('Slope Histogram')
% xlabel('Slope Value')
% ylabel('Counts')
figure (2)
plot(1:s,chain(:,1))
title('Slope Trace')
xlabel('Step')
ylabel('Slope Value')
slope_conf = sqrt(sum((chain(:,1)-t(1)).^2)./s);
% figure (3)
% hist(chain(:,2),100)
% title('Intercept histogram')
% xlabel('Slope Value')
% ylabel('Counts')
figure (4)
inter_conf = sqrt(sum((chain(:,2)-t(2)).^2)./s);
plot(1:s,chain(:,2))
title('Intercept Trace')
xlabel('Step')
ylabel('Intercept Value')

% acceptance = accepts/(accepts+rejects);
% s_binsize = 0.002;
% s_values = 2.45:s_binsize:2.6;
% bin_s = histc(chain(:,1),s_values);
% i_binsize = 0.01;
% i_values = 9.4:i_binsize:10.4;
% bin_i = histc(chain(:,2),i_values);
% post2D = zeros(length(bin_s),length(bin_i));
% for i = 1:length(bin_s)
%     for j = 1:length(bin_i)
% post2D(i,j) = bin_s(i).*bin_i(j);
%     end
% end
% post2D = post2D./sum(sum(post2D*0.002*0.01));
% contour(i_values,s_values,post2D)
% title('2D Posterior')
% xlabel('Intercept')
% ylabel('Slope')

% line = t(1).*X+t(2);
% errorbar(X,Y,sigma,'.')
% hold on
% plot(X,line,'g')
% title('MCMC Best Fit Test')
% xlabel('X values')
% ylabel('Y values')