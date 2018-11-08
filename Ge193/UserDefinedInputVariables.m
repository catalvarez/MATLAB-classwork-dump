function [x,s,S,b,B,H0,AGlen,n,C,C0,m,rho,rhow,g,alpha,nTimeSteps,time,dt,...
    connectivity,nip,Amfcu,blambdau,MFCh,Amfch,blambdah,NumConstants]=UserDefinedInputVariables()

%% nodal variables:
% x     : nodal coordinates
% s     : surface 
% b     : bed
% C     : basal slipperiness
%
%% constants:
% AGlen : rate factor in Glens flow law
% n     : stress exponent in Glens flow law
% m     : stress exponent in basal sliding law
% rho   : ice density
% rhow  : ocean/water density
% g     : gravitational acceleration
% alpha : tilt of the coordinate system (alpha=0 for a vertical z axis)
% nTimeSteps : number of time steps
% time  : model time at the begin of run 
% dt    : time step
%
%% mesh 
% connectivity : FE connectivity
% nip          : number of element integration points
% nod          : number of nodes per element (either 2 or 3, but can (small modifications needed) be up to 5)
%
%% Boundary conditions
% Amfcu    : multi-linear constraint matrix for velocity
% blambdau : the right-hand side of the multi-linear constraint-system for velocity
% Amfch & blambdah : save for thickness
%
%%

H0=100;                  % initial ice thickness
dx=H0; x1=0 ; x2=5000*dx;  % mesh is from x1 to x2 with element size dx
nod=3; nip=3 ;            % Nnodes and integration points per element

[x,connectivity]=CreateMesh1d(x1,x2,dx,nod);

alpha=0.000 ; % alpha is the tilt of the coordinate system
rhow=1000;  % ocean density
%m=3 ; n=3; AGlen=4.5e-24*1e9*24*60*60; rho=917.; g=9.81/1000 ; C0=1/64000;
%m=3 ; n=1; AGlen=4.5e-10; rho=917.; g=9.81/1000 ; C0=1/64000;
m=1 ; n=1; AGlen=1.0e-7; rho=917.; g=9.81/1000 ; C0=1/30;
%alpha=0.005; m=3 ; n=1; AGlen=1.0e-7; rho=917.; g=9.81/1000 ; C0=1/100;





nTimeSteps=10000;  % number of time steps
time=0; % time at beginnig of run
dt=10;          % size of time step, years
sinusoid=0;      % if true then bed is a sinusoid





% define s, b, h and C for each node

%tilted surface
B = -0.005*x+600;   % bedrock geometry
b = B;              % initial ice base
s = b + H0 ;        % initial ice surface
S = 300 ;           % ocean surface
% h = s-b;            % ice thickness
C = x.*0 + C0 ;
% C=C0;           % ice slipperiness



if sinusoid
    lambda=300*1200; db=0.1*H0 ;dC=0;
    
    b=db*sin(2*pi*x/lambda);
    h=s-b;
    C=C0*(1+dC*sin(2*pi*x/lambda));
    save PertParm lambda db H0 C0 dC
end


MFCh=1;  % set to 1 if h constraints enforced using Lagrange mulipl, 2 if using penalty method, 0 if no h constraints in model
% MFC allow for a very general way of defining constrains between nodal values

% A_pq u_q= b_p
Nnodes=numel(x) ; 
% Amfcu=zeros(1,Nnodes);Amfcu(1,1)=1 ; Amfcu(1,end)=-1 ; blambdau=[0]; % u1-un=0  periodic
% Amfch=zeros(1,Nnodes);Amfch(1,1)=1 ; Amfch(1,end)=-1 ; blambdah=[0]; % h1-hn=0  periodic

Amfcu=zeros(1,Nnodes);Amfcu(1,1)=1 ; blambdau=[0];   % u defined at left edge
Amfch=[]; blambdah=[];

%Amfch=zeros(2,length(Nnodes));Amfch(1,1)=1 ; Amfch(2,end)=1 ; blambdah=[H0 ; H0]; % h Dirichlet


%--------------]
NumConstants.dh1tol=1e-4;  % tolerance for change in h using `corrector' in time integration step
NumConstants.du1tol=0.001 ; % tolerance for change in u using `corrector' in time integration step
NumConstants.dtmin=0.001; NumConstants.dtmax=100 ; % minimum and maximum time step
NumConstants.etamax=1e20;
NumConstants.NRtol=1e-8;
NumConstants.Picardtol=1e-7;
NumConstants.MaxNR=50 ; NumConstants.MaxPicard=50;
NumConstants.Uzawatol=1e-7;
NumConstants.GLmeshRefinementLevel=15;
NumConstants.GLmeshRefinementFactor=1.2;
NumConstants.kHe=1e-2;
NumConstants.theta=0.5 ;     % theta=0 is explicit, theta>=0.5 for stability
NumConstants.gamma=1 ;       % 1 for third-order Taylor-Galerkin    
NumConstants.InfoLevel=1; 

%save mesh1d.mat Nnodes nod  x connectivity s b h C nip
%save mparameters alpha m n AGlen rho g nTimeSteps theta dt Amfcu blambdau MFCh Amfch blambdah NumConstants

end