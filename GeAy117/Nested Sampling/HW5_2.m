%Imports data
mat = dlmread('ps5_data2.txt');
P = mat(:,1);
rho  = mat(:,2);
sigma = mat(:,3);

%variables for number of active points and number of samples collected
nActive = 30;
nSamples = 350;

% arrays for storing parameters, areas, and likelihoods
sample_param = zeros(nSamples,3);
sample_logL = zeros(nSamples,1);
sample_logA = zeros(nSamples,1);

% starting area and width
logZ = -9999999999999;
current_logwidth = log(1-exp(-1./nActive));

% drawing each active sample uniformly from their prior spaces
A_K_prime = [1:30] ./30.*3+3;
A_K0 = [1:30] ./30.*30+120;
A_rho0 = [1:30] ./30.*1+3;
A_param = [A_K_prime; A_K0; A_rho0];

for i = 1:nSamples
    % compute likelihoods for active samples
    A_logL = get_log_likelihood(A_param(3,:),A_param(2,:),A_param(1,:),P,rho,sigma);
    % find worst likelihood value
    [log_Lstar,index] = min(A_logL);
    % compute area of worst point
    logA = current_logwidth+log_Lstar;
    % add worst sample to the storage arrays
    sample_param(i,:) = [A_param(1,index) A_param(2,index) A_param(3,index)];
    sample_logL(i) = log_Lstar;
    sample_logA(i) = logA;
    % update total area
    logZ = update_logZ(logZ,logA);
    % draw a new sample
    [new_param,logL_new] = get_new_sample(A_param,A_logL,P,rho,sigma);
    % overwrite worst point with the new sample
    A_param(:,index) = new_param;
    % update the width
    current_logwidth = current_logwidth - 1./nActive;
end

Z = exp(logZ);
A = exp(sample_logA);
W = A./Z;

Ave_K_prime = sum(sample_param(:,1).*W);
Ave_K0 = sum(sample_param(:,2).*W);
Ave_rho0 = sum(sample_param(:,3).*W);

figure (1)
line = Ave_rho0.*(1+Ave_K_prime.*P./Ave_K0).^(1./Ave_K_prime);
errorbar(P,rho,sigma,'.')
hold on
plot(P,line,'g')
title('Nested Sampling')
xlabel('P (GPa)')
ylabel('density (g/cm^3)')

figure (2)
subplot(2,2,1)
plot(sample_param(:,1),sample_param(:,3),'.')
title('Density vs dK/dP')
xlabel('K prime')
ylabel('rho')
subplot(2,2,2)
plot(sample_param(:,2),sample_param(:,3),'.')
title('Density vs Bulk Modulus')
xlabel('K_0')
ylabel('rho')
subplot(2,2,3)
plot(sample_param(:,2),sample_param(:,1),'.')
title('dK/dP vs Bulk Modulus')
xlabel('K_0')
ylabel('K prime')
