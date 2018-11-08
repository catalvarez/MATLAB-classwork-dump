function [u,dudx,xint,etaInt,kv,rh,lmfc]=SSHELF1d(nels,nod,nodes,gcoord,connectivity,nip,...
    s,b,...
    iteration_max,iteration_min,tol,etamax,etaInt,u0,lmfc,...
    alpha,g,rho,AGlen,beta,n,m,...
    Amfc,blambda,Neumann,dudxb,dudxa)


% solves the non-linear 1d shallow ice stream equation
% using the method of finite elements
%
%  nf
rhogcos=rho*g*cos(alpha); rhogsin=rho*g*sin(alpha);

ndim=1;
dof=1;  nodes_fixed=[]; % the steering vector holds the freedom numbers for each element g(nod*dof,nels)

[~,g_g]=create_steering_vector(connectivity,nodes,nodes_fixed,dof);

neq=max(g_g(:));

h=s-b;

[points,weights]=sample(ndim,nip,nod);  % get local coordinates and weights
% if using different elements, this
% needs to be within element loop


% non-linear iteration loop

diff=1e10;
if n==1 && m==1 ; iteration_max=1 ; iteration_min=1 ; end
iteration=0;
while (diff > tol  && iteration < iteration_max )  || iteration < iteration_min
    iteration=iteration+1;
    rh=zeros(neq,1) ; kv=zeros(neq,neq)  ;
    % element loop
    for Iele=1:nels
        % gather local quantities from global arrays
        
        coord=gcoord(connectivity(Iele,2:end));
        h_l=h(connectivity(Iele,2:end));
        s_l=s(connectivity(Iele,2:end));
        beta_l=beta(connectivity(Iele,2:end));
        
        g_l=g_g(:,Iele); % contains the freedom number associated with the element
        u0_l=u0(connectivity(Iele,2:end));
        
        km=zeros(nod,nod); rhs=zeros(nod,1);
        for Iint=1:nip                           % loop over integration points
            
            fun=shape_fun(Iint,ndim,nod,points);  % nod x 1  : values of form functions at integration points
            der=shape_der(Iint,ndim,nod,points);  % 1 x nod  : dNdXi (derivatives of formfunctions with respect to local coordinates at inegration poitns)
            J=der*coord;                           % Jacobian
            %iJ=inv(J);
            deriv=der/J;                          % derivatives of form functions with respecty to global coordinates
            D=-4*etaInt(Iele,Iint)*(h_l'* fun)           ;
            lhs1=D* (deriv' * deriv)*det(J)*weights(Iint);  % 1st term on left-hand side (Stiffness)
            
            %Ceff=C_l.^(1/m).*u0_l.^(1-1/m);
            %lhs2=-   fun*fun'*det(J)*weights(Iint)/(Ceff'*fun);   % 2nd term on left-hand side
            
            
            betaeff2=beta_l.^(2/m).*u0_l.^(1/m-1);
           
            lhs2=-fun*fun'*det(J)*weights(Iint)*(betaeff2'*fun);   % 
            
            km=km+lhs1+lhs2;
            
            rhs1=rhogcos* (h_l'*fun) * s_l'*deriv'*fun*det(J)*weights(Iint);
            rhs2=-rhogsin* (h_l'*fun) *fun*det(J)*weights(Iint);
            rhs=rhs+rhs1+rhs2;
            
        end % integration points
        
        
        % assemble global matrix
        for i1=1:length(g_l)
            for i2=1:length(g_l)
                kv(g_l(i1),g_l(i2))=kv(g_l(i1),g_l(i2))+km(i1,i2);
            end
            rh(g_l(i1))=rh(g_l(i1))+rhs(i1);
        end
        
       
        
    end  % element loop
    
    
    % [-  Boundary conditions
    % ( Neumann assuming that BC are only for element 1 and nels)
    if Neumann
        if ~isnan(dudxb)            
            Iele=nels;
            rhs=-4*etaInt(Iele,end)*h(g_l(end))*dudxb;
            rh(g_l(end))=rh(g_l(end))+rhs;
        end
        
        if ~isnan(dudxa)
            % Neumann at left-hand side
            Iele=1;
            rhs=-4*etaInt(Iele,1)*h(g_l(1))*dudxa;
            rh(g_l(1))=rh(g_l(1))+rhs;
            
        end
    end
    %[- MultiFreedom constraints (MFCs)
    % deal with MFCs using Lagrange multipliers
    
    
    [u,lmfc]=solveKApeSymmetric(kv,Amfc,rh,blambda,lmfc);

    
    diff=abs(max(u-u0)/max(u)); % vel increment/ max vel
    %disp([' |max(u-u0)/max(u)| = ',num2str(diff)])
    
    %
    % get strain rates at integration points
    %dudx=dNdXi dXidx up
    %
    dudx=zeros(nels,nip); xint=zeros(nels,nip); 
    for Iele=1:nels                            % loop over elements
        g_l=g_g(:,Iele);
       
        u_l=u(connectivity(Iele,2:end)); % velocity at nodal point of element
        coord=gcoord(connectivity(Iele,2:end));  % coordinates of nodal points
        for Iint=1:nip                           % loop over integration points
            
            fun=shape_fun(Iint,ndim,nod,points);  % nod x 1, values of form functions at integration points
            der=shape_der(Iint,ndim,nod,points);   % dNdXi
            J=der*coord;                           % Jacobian
            %iJ=inv(J);
            deriv=der/J;                          % dNdx
            
            
            dudx(Iele,Iint)=deriv*u_l;          % dudx at integration point Iint of element Iele
            xint(Iele,Iint)=fun'*coord ;
            etaInt(Iele,Iint)=min([0.5*AGlen^(-1/n)*(abs(dudx(Iele,Iint)))^((1-n)/n),etamax]);
            
            
        end
        
    end
    
    %figure(10) ; hold on ;plot(xint(:),etaInt(:),'o'); title(' etaInt')
end   % non-linear iteration loop



