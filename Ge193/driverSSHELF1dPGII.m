clear all
%%
% A simple wrapper around repeated calls to solve SSA and the thickness evolution equation
% Solves the SSTREAM/SSA equation in 1d using FE
% Solves the thickenss evolution using third-order Taylor-Galerkin
%

%% Get geometry, rheology, etc as defined by user
 [x,s,S,b,B,H0,AGlen,n,C,C0,m,rho,rhow,g,alpha,nTimeSteps,time,dt,...
    connectivity,nip,Amfcu,blambdau,MFCh,Amfch,blambdah,NumConstants]=UserDefinedInputVariables();

%% Initialize/define some variables
Nnodes=numel(x) ; [Nele,nod]=size(connectivity);
h=s-b;   % ice thickness
lambdau=blambdau*0;
u=s*0; dudt=s*0; % needed for prognostic step
lmfc=blambdau*0;
gfele=zeros(Nele,1)+1;
S=S+s*0;
H = S-B;

plotfig=1;   % set to one if figures to be plotted
calc_transfunc=0;  % calculate analytical trans func for comparision with num. results


%% If interested in comparing with analytical transfer functions keep this part
% otherwise just delete
% if n~=1 ; calc_transfunc=0 ; H0=1; end
% if calc_transfunc
% 	load PertParm lambda db H0 C0 dC
% 	uAmplAnalytical=zeros(max([1,nTimeSteps]),1); uAmplNumerical=uAmplAnalytical;
% 	uPhaseAnalytical=uAmplAnalytical; uPhaseNumerical=uAmplAnalytical;
% 	sAmplAnalytical=zeros(max([1,nTimeSteps]),1); sAmplNumerical=sAmplAnalytical;
% 	sPhaseAnalytical=sAmplAnalytical; sPhaseNumerical=sAmplAnalytical;
% 	
% end
%% Transient calculation
% uses staggered/semi-implicit approach
% 1) calculate velocity (u) as a function of geometry, etc   (diagnostic step)
% 2) calculate thickenss h at the end of a time step for given velocity, geometry, etc (prognostic step)
%



for Itime=1:max([1,nTimeSteps])
        I = rho .* h < rhow .* H ; 
        s(I) = S(I) +(1-rho/rhow).*h(I); % new floating ice surface
        s(~I) = B(~I) + h(~I) ;             % new grounded ice surface
      C(~I) = C0;
        C(I) = 10^10;
    
	ulast=u;

    % diagnostic step
	 [u,lambdau,kv,rh]=SSTREAM1d(s,h,u,x,connectivity,nip,gfele,AGlen,C,Amfcu,blambdau,lambdau,n,m,alpha,rho,rhow,g,NumConstants);
	 
	
	dudtm1=dudt;  % needed for prognostic step
	dudt=(u-ulast)/dt;

%% some plotting, change this as needed	
	if plotfig
        if ~mod(Itime,(nTimeSteps/1000))
% 		figure(3); hold on ;plot(x/H0,u) ; title(' u ') ; xlabel(' x/H0 ')
% 		figure(2); hold on; plot(x/H0,h,'b') ; title('Thickness') ; xlabel(' x/H0 ')
		figure(1);
        hold on;
        plot(x/H0,s,'r') ;
        plot(x/H0,b,'c') ;
        plot(x/H0,S,'b') ;
        plot(x/H0,B,'k') ;
        legend('ice surface','ice base','ocean','bedrock')
        title('Surf') ;
        xlabel(' x/H0 ') ;
        end
%         ylim([-800 1500]);
        
% 		figure(4); hold on ;plot(x/H0,b,'c')  ; title('Bed') ; xlabel(' x/H0 ')
    end
    
    %% comparision with analytical transfer functions (just delete this if of no interest)
	if calc_transfunc
		% all variables dimensional
		% only makes sense for n=1, but m can have any value
		
		eta=1/2/AGlen; % viscosity for n=1
		
		% horizontal velocity
		tub=SIS_Tub_t_3d_m(2*pi/lambda,eps,time,alpha,H0,eta,C0,rho,g,m);
		tuc=SIS_Tuc_t_3d_m(2*pi/lambda,eps,time,alpha,H0,eta,C0,rho,g,m);
		tu=db*tub+dC*tuc;
		
		[uphase, uampl]=harm(u-mean(u),2*pi/lambda,x);
		auphase=-atan2(imag(tu),real(tu))*180/pi;
		
		% surface topography
		tsb=SIS_Tsb_t_3d_m(2*pi/lambda,eps,time,alpha,H0,eta,C0,rho,g,m);
		tsc=SIS_Tsb_t_3d_m(2*pi/lambda,eps,time,alpha,H0,eta,C0,rho,g,m);
		ts=db*tsb+dC*tuc;
		[sphase, sampl]=harm(s-mean(s),2*pi/lambda,x);
		asphase=-atan2(imag(ts),real(ts))*180/pi;
		
		
		%disp(['  Nondimensional Slipperiness : ',num2str(C0*eta*2/H0)])
		disp(['                         time : ',num2str(time)])
		
		disp(['  u analytical/numerical ampl : ',num2str(abs(tu),8),'/',num2str(uampl,8)])
		disp([' u analytical/numerical phase : ',num2str(auphase),'/',num2str(uphase)])
		
		
		disp(['  s analytical/numerical ampl : ',num2str(abs(ts),8),'/',num2str(sampl,8)])
		disp([' s analytical/numerical phase : ',num2str(asphase),'/',num2str(sphase)])
		disp('----------------------------------------------------------------------------')
		
		uAmplAnalytical(Itime)=abs(tu); uAmplNumerical(Itime)=uampl;
		uPhaseAnalytical(Itime)=auphase; uPhaseNumerical(Itime)=uphase;
		sAmplAnalytical(Itime)=abs(ts); sAmplNumerical(Itime)=sampl;
		sPhaseAnalytical(Itime)=asphase; sPhaseNumerical(Itime)=sphase;
		
		
    end
    
	%% prognostic step
	if ~nTimeSteps == 0
		
		h0=h ; u0=u ; u1=u ; a0=x.*0+0.3 ; a1=x.*0+0.3 ; %a0=u1*0 ; a1=a0 ;
        % get an (explicit) estimate for u at the end of the time step 
        if Itime==1
            u1=u0;
        elseif Itime==2
            u1=u0+dt*dudt;
        else
            u1=u0+dt*( 3*dudt-dudtm1)/2  ; % Adams-Bashforth method
        end

         % calculate h at the end of the time step
         h1=nexthTG3(h0,u0,u1,dudt,a0,a1,dt,NumConstants.theta,NumConstants.gamma,x,connectivity,nip,MFCh,Amfch,blambdah);

    end
    
	 % update variables
        
        b = s-h ;
        h=h1;
		time=time+dt;	
       
     % find x value of grounding line
     % Xgl[Itime] = 
    
	end
	
	
	
	

%% comparision with transfer functions 
if calc_transfunc
	figure(10);
	plot([0:nTimeSteps-1]*dt,uAmplAnalytical,'r-+') ; hold on
	plot([0:nTimeSteps-1]*dt,uAmplNumerical,'g-o') ; xlabel('time'); legend('Analytical','Numerical'); title('u Amplitudes')
	figure(11);
	plot([0:nTimeSteps-1]*dt,uPhaseAnalytical,'r-+') ; hold on
	plot([0:nTimeSteps-1]*dt,uPhaseNumerical,'g-o') ; xlabel('time'); legend('Analytical','Numerical'); title('u Phase')
	
	figure(12);
	plot([0:nTimeSteps-1]*dt,sAmplAnalytical,'r-+') ; hold on
	plot([0:nTimeSteps-1]*dt,sAmplNumerical,'g-o') ; xlabel('time'); legend('Analytical','Numerical'); title('s Amplitudes')
	figure(13);
	plot([0:nTimeSteps-1]*dt,sPhaseAnalytical,'r-+') ; hold on
	plot([0:nTimeSteps-1]*dt,sPhaseNumerical,'g-o') ; xlabel('time'); legend('Analytical','Numerical'); title('s Phase')
end

% disp(lambda)
disp(min(b))
disp(max(b))
% disp(sampl)