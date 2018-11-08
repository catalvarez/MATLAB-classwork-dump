function Daisy_out = daisyworld(alpha_w_i,alpha_b_i,T_i,ntime,alpha_k_i)
%DAISYWORLD is a model based on Watson & Lovelock (1983).
% alpha_w_i: Fraction of land covered by white daisies
% alpha_b_i: Fraction of land covered by black daisies
% T_i:       Initial temperature of daisy world in K (default = 300)
% ntime:     Number of timesteps in main loop (default = 100)
% alpha_k_i: Fraction of land covered by kudzu
% 
%daisyworld.m gka 29 March 2008%
%modified aj March 31 2009%
%modified kb March 25 2010%
%modified ks Dec   31 2013%

%PARAMETERS%
T_opt_w = 295; %Optimal growth temperature of white daisies, in K 
T_opt_b = 295; %Optimal growth temperature of black daisies, in K 
T_opt_k = 308; %Optimal growth temperature of kudzu, in K 
k = 17.5^(-2); %Growth rate bracketed between 35 K of T_opt
albedo_g = 0.5;  %Albedo of bare ground
albedo_w = 0.75; %Albedo of white daisies
albedo_b = 0.25; %Albedo of black daisies
albedo_k = 0.3;  %Albedo of kudzu
S = 1368/4; %Solar radiation in W/m2 received by a planet located 1 au from the Sun
L = 3; %Luminosity (L is set to 3 so that the radiative equilibrium temp of the planet mimics the Earth's)
sigma = 5.670373*10^-8;% Stefan Boltzman constant W/m2/K4
q = 0.2*S*L/sigma; %Heat transfer coefficient to ensure thermal equilibrium (set so that q < 0.2*SL/sigma)
gamma = 0.3; %Death rate (default = 0.3)

%INITIALIZATION%
if ~exist('alpha_w_i','var') || isempty(alpha_w_i); alpha_w_i = 0; end 
if ~exist('alpha_b_i','var') || isempty(alpha_b_i); alpha_b_i = 0; end 
if ~exist('T_i','var') || isempty(T_i); T_i = 300; end 
if ~exist('ntime','var') || isempty(ntime); ntime = 50; end; ntime = ntime+1;
alpha_w = zeros(1,ntime); %Storage vector for alpha_w
alpha_b = zeros(1,ntime); %Storage vector for alpha_b
alpha_g = zeros(1,ntime); %Storage vector for alpha_b
T = zeros(1,ntime); %Storage vector for Temperature
if alpha_w_i+alpha_b_i>1; error('You have daisies covering more than the available land area.'); end

% Model %
switch nargin
 case 5 % Keeling Catastrophe Model
     alpha_k = zeros(1,ntime); %Storage vector for alpha_k
     L = 0.75*L; % Reduce luminosity to 75% of that of the Earth
     q = 0.2*S*L/sigma;
     
 % Initialize with black daises until they propagate the barren planet
  for itime = 1:50
	alpha_g = 1 - alpha_b_i; %Fraction of land that is bare ground
        if alpha_g<0; error('You have daisies covering more than the available land area.'); end
	A = alpha_b_i*albedo_b + albedo_g*(alpha_g); %mean planetary albedo
	alpha_b(itime) = alpha_b_i;
	T(itime) = T_i;
    
	T_i = (S*L*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
	T_b = (q*(A - albedo_b) + T_i^4)^(1/4); %compute temperature of patch of black daisies
	alpha_b_i = alpha_b_i + alpha_b_i*(alpha_g*betafn(T_b,T_opt_b,k) - gamma);
  end 
  figure;
  subplot(2,3,1); plot(0:49,T(1:50),'r-'); xlim([0 50]); ylabel('Temp (K)'); title('Sequence 1')
  subplot(2,3,4); plot(0:49,alpha_b(1:50),'b-'); xlim([0 50]);
  	ylim([0 1]); ylabel('fractional cover'); legend('black')
    
 % Introduce kudzu until black daisies die off
  for itime = 51:151
    alpha_g = 1 - alpha_b_i - alpha_k_i; %Fraction of land that is bare ground
        if alpha_g<0; error('You have daisies covering more than the available land area.'); end
	A = alpha_b_i*albedo_b + alpha_k_i*albedo_k + albedo_g*(alpha_g); %mean planetary albedo
	alpha_b(itime) = alpha_b_i;
    alpha_k(itime) = alpha_k_i;
	T(itime) = T_i;
      
	T_i = (S*L*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
	T_b = (q*(A - albedo_b) + T_i^4)^(1/4); %compute temperature of patch of black daisies
    T_k = (q*(A - albedo_k) + T_i^4)^(1/4); %compute temperature of patch of kudzu
      
    alpha_b_i = alpha_b_i + alpha_b_i*(alpha_g*betafn(T_b,T_opt_b,k) - gamma);
    alpha_k_i = alpha_k_i + alpha_k_i*(alpha_g*betafn(T_k,T_opt_k,k) - gamma);  
  end
    subplot(2,3,2); plot(0:100,T(51:151),'r-'); xlim([0 100]); ylabel('Temp (K)'); title('Sequence 2')
    subplot(2,3,5); plot(0:100,alpha_b(51:151),'b-',0:100,alpha_k(51:151),'m-') 
        xlim([0 100]); ylim([0 1]); legend('black','kudzu')
    
 % Introduce white daisies and witness the annihilation
  for itime = 152:751
    alpha_g = 1 - alpha_w_i - alpha_k_i; %Fraction of land that is bare ground
        if alpha_g<0; error('You have daisies covering more than the available land area.'); end
	A = alpha_w_i*albedo_w + alpha_k_i*albedo_k + albedo_g*(alpha_g); %mean planetary albedo
	%alpha_b(itime) = alpha_b_i;
    alpha_k(itime) = alpha_k_i;
	alpha_w(itime) = alpha_w_i;
	T(itime) = T_i;
    
	T_i = (S*L*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
	%T_b = (q*(A - albedo_b) + T_i^4)^(1/4); %compute temperature of patch of black daisies
    T_k = (q*(A - albedo_k) + T_i^4)^(1/4); %compute temperature of patch of kudzu
    T_w = (q*(A - albedo_w) + T_i^4)^(1/4); %compute temperature of patch of white daisies
      
    %alpha_b_i = alpha_b_i + alpha_b_i*(alpha_g*betafn(T_b,T_opt_b,k) - gamma);
    alpha_k_i = alpha_k_i + alpha_k_i*(alpha_g*betafn(T_k,T_opt_k,k) - gamma); 
	alpha_w_i = alpha_w_i + alpha_w_i*(alpha_g*betafn(T_w,T_opt_w,k) - gamma);
  end
    subplot(2,3,3); plot(0:599,T(152:751),'r-'); xlim([0 600]); ylabel('Temp (K)'); title('Sequence 3')
    subplot(2,3,6); plot(0:599,alpha_k(152:751),'m-',0:599,alpha_w(152:751),'c-') 
        xlim([0 600]); ylim([0 1]); legend('kudzu','white')
    set(gcf,'color','w');
    
 otherwise
  for itime = 1:ntime
	alpha_g(itime) = 1 - alpha_b_i - alpha_w_i; %Fraction of land that is bare ground
        if alpha_g(itime)<0; error('You have daisies covering more than the available land area.'); end
    A = alpha_w_i*albedo_w + alpha_b_i*albedo_b + albedo_g*(alpha_g(itime)); %mean planetary albedo
	if ntime==10^5+1; S=1368/4*(1+0.3*sin(2*pi*itime/10^4)); end %insolation for problem 3
    
	%store variables for plotting%
    alpha_w(itime) = alpha_w_i; 
	alpha_b(itime) = alpha_b_i;
	T(itime) = T_i;
    
	T_i = (S*L*(1-A)/sigma)^(1/4); %compute mean planetary temperature in radiative equilibrium
	T_w = (q*(A - albedo_w) + T_i^4)^(1/4); %compute temperature of patch of white daisies
	T_b = (q*(A - albedo_b) + T_i^4)^(1/4); %compute temperature of patch of black daisies
	T_g = (q*(A - albedo_g) + T_i^4)^(1/4); %compute temperature of bare ground
	
    % For fun: Introduce evolutionary changes to model
    %if mod(itime,10^2)==1     
    %   T_opt_b = T_opt_b + 0.1*(T - T_opt_b);
	%end

	alpha_w_i = alpha_w_i + alpha_w_i*(alpha_g(itime)*betafn(T_w,T_opt_w,k) - gamma);
	alpha_b_i = alpha_b_i + alpha_b_i*(alpha_g(itime)*betafn(T_b,T_opt_b,k) - gamma);
  end
  
  %DATA PLOTTING%
    figure%(1);
    subplot(2,1,1); plot(0:ntime-1,T,'r-')
    ylabel('Temp (K)')
    subplot(2,1,2); plot(0:ntime-1,alpha_b,'b-',0:ntime-1,alpha_w,'c-') 
    ylim([0 1])
    ylabel('fractional cover')
    legend('black', 'white')
    set(gcf,'color','w');
end

% Save Variables to Workspace % 
V = who;
for i = 1:length(V)
    Daisy_out.(V{i}) = eval(V{i});
end

function beta = betafn(temp, T_opt, k)
  if abs(T_opt - temp) < k^(-1/2)
    beta = 1 - k*(T_opt - temp)^2;
  else
    beta = 0;
  end
return