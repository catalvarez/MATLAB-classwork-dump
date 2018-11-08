function [u,uitr,dudx,xint,etaInt,kv,rh,sol,indu]=ShallowStreamShelf1dBC_PGII(nels,nod,nodes,gcoord,connectivity,nip,...
        s,h,gfele,...
        iteration_max,iteration_min,tol,etamax,etaInt,ulast,...
        alpha,g,rho,rhow,AGlen,beta,n,m,...
        MFC,Amfc,blambda)
    
    persistent SF DR eqnr g_g neq points weights
        
    % right hand-side modified so that Weertman BC are the natural boundary conditions
    % written and tested on 26 August, 2009
    
    
    % I need a realistic upper value for beta2
    usmall=0.01;  % velocity that is small compared to most expected velocities
    beta2max=zeros(nod,1)+max(beta(:))^(2/m)*usmall^(1/m-1);
        
    uitr=ulast;
    
    % solves the non-linear 1d shallow ice stream equation
    % and shallow shelf equation
    % using the method of finite elements
    %
    %  nf
    rhogcos=rho*g*cos(alpha); rhogsin=rho*g*sin(alpha);
    
    ndim=1;
    dof=1;  nodes_fixed=[]; % the steering vector holds the freedom numbers for
    % each element g(nod*dof,nels)
    % i.e. the labelling of variables in the
    % system of equations
    
    if isempty(SF)
        [eqnr,g_g]=create_steering_vector(connectivity,nodes,nodes_fixed,dof);
        neq=max(g_g(:));
        [points,weights]=sample(ndim,nip,nod);  % get local coordinates and weights
        % if using different elements, this
        % needs to be within element loop
        SF=zeros(nip,nod);DR=SF;
        for Iint=1:nip
            SF(Iint,1:nod)=shape_fun(Iint,ndim,nod,points);  % nod x 1  : values of form functions at integration points
            DR(Iint,1:nod)=shape_der(Iint,ndim,nod,points);  % 1 x nod  : dNdXi (derivatives of formfunctions with respect to local coordinates
        end
    end
    % non-linear iteration loop
    
    diff=1e10;
    if n==1 && m==1 ; iteration_max=1 ; iteration_min=1 ; end
    iteration=0;
    while (diff > tol  && iteration < iteration_max )  || iteration < iteration_min
        iteration=iteration+1;
        %disp([' SSS iteration # : ',num2str(iteration)])
        rh=zeros(neq,1) ; kv=zeros(neq,neq)  ;
        
        % element loop
        for Iele=1:nels
            % gather local quantities from global arrays
            
            
            coord=gcoord(connectivity(Iele,2:end));
            h_l=h(connectivity(Iele,2:end));
            s_l=s(connectivity(Iele,2:end));
            beta_l=beta(connectivity(Iele,2:end));
            
            g_l=g_g(:,Iele); % contains the freedom number associated with the element
            uitr_l=uitr(connectivity(Iele,2:end));
            
            km=zeros(nod,nod); rhs=zeros(nod,1);
            
            % check if the element has nodes downstream of the gl
            
            gf_l=gfele(Iele);
            %if gf_l ==0 ; disp(['element floating ',num2str(Iele)]); end
            
            for Iint=1:nip                           % loop over integration points
                
                %fun=shape_fun(Iint,ndim,nod,points);  % nod x 1  : values of form functions at integration points
                %der=shape_der(Iint,ndim,nod,points);  % 1 x nod  : dNdXi (derivatives of formfunctions with respect to local coordinates at inegration poitns)
                
                fun=SF(Iint,:)'; der=DR(Iint,:);
                
                J=der*coord;
                detJ=det(J);
                
                iJ=inv(J);
                deriv=iJ*der;                          % derivatives of form functions with respecty to global coordinates
                D=-4*etaInt(Iele,Iint)*(h_l'* fun)           ;
                lhs1=D* (deriv' * deriv)*det(J)*weights(Iint);  % 1st term on left-hand side (Stiffness)
                
                %Ceff=C_l.^(1/m).*uitr_l.^(1-1/m);
                %lhs2=-   fun*fun'*detJ*weights(Iint)/(Ceff'*fun);   % 2nd term on left-hand side
                
                betaeff2=beta_l.^(2/m).*abs(uitr_l).^(1/m-1);
                if m~=1
                    betaeff2=min([betaeff2' ; beta2max']); betaeff2=betaeff2';
                end
                
                lhs2=-gf_l*(fun*fun')*detJ*weights(Iint)*(betaeff2'*fun);   %
                
                km=km+lhs1+lhs2;
                
                
                rhs1=rhogcos* (h_l'*fun) * (s_l'*deriv'-(1-rho/rhow)* h_l'*deriv') *fun; 
                %rhs1=rhogcos* (h_l'*fun) * s_l'*deriv'*fun;
                
                rhs2=-rhogsin* (h_l'*fun) *fun;
                
                rhs3=-0.5*rhogcos*(1-rho/rhow)* (h_l'*fun)^2*deriv';
                rhs=rhs+(rhs1+rhs2+rhs3)*detJ*weights(Iint);
                
            end % integration points
            
            
            % assemble global matrix
            for i1=1:length(g_l)
                for i2=1:length(g_l)
                    kv(g_l(i1),g_l(i2))=kv(g_l(i1),g_l(i2))+km(i1,i2);
                end
                rh(g_l(i1))=rh(g_l(i1))+rhs(i1);
            end
            
            
            
        end  % element loop
        

        if MFC==0  % no mfc constrains
            Amfc=[] ;blambda=[];
        end
        
        
        if MFC==3   % take unknowns from matrix
            Amfc=[] ; blambda=[];
            kv2=kv(2:end,2:end) ; rh2=rh(2:end) ; sol=real(sparse(kv2)\rh2);
            u=[eps ; sol]; indu=[]; lmfc=[];
            
        else
            %disp(' lagrange ')
            [nA,mA]=size(Amfc);
            kv=[kv Amfc' ; Amfc zeros(nA,nA)]; rh=[rh ; blambda] ;
            sol=real(sparse(kv)\rh);
            indu=find(eqnr(:,2)==1) ;
            u(eqnr(indu,1))=sol(indu); u=u(:);
            lmfc=sol(end-nA+1:end);
        end
        
        
        diff=abs(max(u-uitr)/mean(abs(u))); % vel increment/ max vel
        %disp([' Iteration : ',num2str(iteration),' |max(u-uitr)/max(u)| = ',num2str(diff)])
        
        uitr=u;
        
        %
        % get strain rates at integration points and calculate new effective viscosity
        %dudx=dNdXi dXidx up
        %
        
        dudx=zeros(nels,nip); xint=zeros(nels,nip); dLdc_l=zeros(nod,nod);
        for Iele=1:nels                            % loop over elements
            g_l=g_g(:,Iele);
            
            
            u_l=u(connectivity(Iele,2:end)); % velocity at nodal point of element
            %u_l=uitr(connectivity(Iele,2:end)); % velocity at nodal point of element
            coord=gcoord(connectivity(Iele,2:end));  % coordinates of nodal points
            for Iint=1:nip                           % loop over integration points
                
                %fun=shape_fun(Iint,ndim,nod,points);  % nod x 1, values of form functions at integration points
                %der=shape_der(Iint,ndim,nod,points);   % dNdXi
                fun=SF(Iint,:)'; der=DR(Iint,:);
                J=der*coord;                           % Jacobian
                iJ=inv(J);
                deriv=iJ*der;                          % dNdx
                
                
                dudx(Iele,Iint)=deriv*u_l;          % dudx at integration point Iint of element Iele
                xint(Iele,Iint)=fun'*coord ;
                etaInt(Iele,Iint)=min([0.5*AGlen^(-1/n)*(abs(dudx(Iele,Iint)))^((1-n)/n),etamax]);
                
                
            end
            
            
        end
    
        
        
    end   %
    
    
    
    