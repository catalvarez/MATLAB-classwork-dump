function [R,T,F]=rh1NR(s,h,u,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g)
    
    % calculates the FE right-hand side in a vectorized form
    
    Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
    ndim=1; dof=1; neq=dof*Nnodes;
    
    
    
    
    
    hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
    snod=reshape(s(connectivity,1),Nele,nod);
    unod=reshape(u(connectivity,1),Nele,nod);
    Cnod=reshape(C(connectivity,1),Nele,nod);
    
    
    
    rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
    
    [points,weights]=sample(ndim,nip,nod);
    
    
    Tx=zeros(Nele,nod); Fx=zeros(Nele,nod);
    
    
    [etaInt]=calcStrainRatesEtaInt1D(u,coordinates,connectivity,nip,AGlen,n);
    
    for Iint=1:nip
        
        
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        [Deriv,detJ]=derivVector1d(coordinates,connectivity,nip,Iint);
        
        
        
        % Deriv : Nele x dof x nod
        %  detJ : Nele
        
        % values at integration this point
        
        % gfint=gfnod*fun;
        % gfint=gfele
        
        
        hint=hnod*fun;
        uint=unod*fun;
        Cint=Cnod*fun;
        
        beta2int = calcBeta2in1Dint(uint,Cint,m); % it would be fairly easy to implement GF.int here
        beta2int=beta2int.*gfele;  % all or none integration points of a element grounded
        
        
        etaint=etaInt(:,Iint) ;  % I could consider calculating this here
        
        
        % derivatives at this integration point for all elements
        dsdx=zeros(Nele,1); dhdx=zeros(Nele,1); dudx=zeros(Nele,1);
        
        for Inod=1:nod
            
            dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
            dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
            dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
            
            
        end
        
        
        detJw=detJ*weights(Iint);
        
        
        for Inod=1:nod
            
            
            %t1=(rhog*hint.*sa-beta2int.*uint-rhog*hint.*ca.*(dsdx-(1-rho/rhow).*dhdx)).*fun(Inod);
            %t2=(0.5*ca*rhog.*(1-rho/rhow).*hint.^2-4*hint.*etaint.*dudx).*Deriv(:,1,Inod);
            
            t1=(rhog*hint.*sa-rhog*hint.*ca.*(dsdx-(1-rho/rhow).*dhdx)).*fun(Inod)+(0.5*ca*rhog.*(1-rho/rhow).*hint.^2).*Deriv(:,1,Inod);  % F
            t2=-beta2int.*uint.*fun(Inod)-4*hint.*etaint.*dudx.*Deriv(:,1,Inod);                                            % T
            
            
            %b1(:,Inod)=b1(:,Inod)+(t1+t2).*detJw;
            
            Tx(:,Inod)=Tx(:,Inod)-t2.*detJw;  % internal nodal forces
            Fx(:,Inod)=Fx(:,Inod)+t1.*detJw;  % external nodal forces
            
            
            
        end
    end
    
    
    
    
    % assemble right-hand side
    
    T=sparse(neq,1); F=sparse(neq,1);
    for Inod=1:nod
        
        T=T+sparse(connectivity(:,Inod),ones(Nele,1),Tx(:,Inod),neq,1);
        F=F+sparse(connectivity(:,Inod),ones(Nele,1),Fx(:,Inod),neq,1);
        
        
    end
    
    R=T-F; % note that I'm solving K x = -R
    
end



