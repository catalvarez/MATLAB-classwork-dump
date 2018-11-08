
function [K,R,T,F]=kvrh1NR(s,h,u,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g)
	
	% calculates the FE matrix and right-hand side in a vectorized form
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity);
	ndim=1; dof=1; neq=dof*Nnodes;
	
		
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	snod=reshape(s(connectivity,1),Nele,nod);
	unod=reshape(u(connectivity,1),Nele,nod);
	Cnod=reshape(C(connectivity,1),Nele,nod);
	
	
	rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
	
	[points,weights]=sample(ndim,nip,nod);
		
	
	K=sparse(neq,neq);
	d1d1=zeros(Nele,nod,nod);
	Tx=zeros(Nele,nod); Fx=zeros(Nele,nod);
	
    
	[etaInt,xint,exx,DetaDuInt]=calcStrainRatesEtaInt1D(u,coordinates,connectivity,nip,AGlen,n);
	
	
	%figure ; plot(xint(:),DetaDuInt(:),'ro') ; title('DetaDuInt')
	%pause(1)

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
		
		
        [beta2int,Dbeta2Duuint] = calcBeta2in1Dint(uint,Cint,m); % it would be fairly easy to implement GF.int here
		Dbeta2Duuint=Dbeta2Duuint.*gfele;
		beta2int=beta2int.*gfele;  % all or none integration points of a element grounded
		
		etaint=etaInt(:,Iint) ;  % I could consider calculating this here
		detaduint=DetaDuInt(:,Iint);
		
		% derivatives at this integration point for all elements
		dsdx=zeros(Nele,1); dhdx=zeros(Nele,1); dudx=zeros(Nele,1);
		
		for Inod=1:nod
			
			dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
			dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
			dudx=dudx+Deriv(:,1,Inod).*unod(:,Inod);
			
			
		end
				
		
		detJw=detJ*weights(Iint);
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+...
					4*hint.*detaduint.*dudx.*Deriv(:,1,Inod).*Deriv(:,1,Jnod).*detJw+...
					4*etaint.*hint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod).*detJw+...
					beta2int.*fun(Jnod).*fun(Inod).*detJw+...
					Dbeta2Duuint.*fun(Inod).*fun(Jnod).*detJw;
				
			end
			
			%t1=-rhog*hint.*((dsdx-(1-rho/rhow).*dhdx)*ca-sa).*fun(Inod);
			%t2=0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,1,Inod); % correction on 11 march 2010
% 			
% 			t1=(rhog*hint.*sa-beta2int.*uint-rhog*hint.*ca.*(dsdx-(1-rho/rhow).*dhdx)).*fun(Inod);
% 			t2=(0.5*ca*rhog.*(1-rho/rhow).*hint.^2-4*hint.*etaint.*dudx).*Deriv(:,1,Inod);
% 			
% 			b1(:,Inod)=b1(:,Inod)+(t1+t2).*detJw;
			
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
    

	
% 	figure(2002) ; plot(rh) ; title('rh')
% 	figure(2003) ; plot(rh(2:end-1)) ; title('rh(2:end-1)')
% 	pause(1)
	% assemble matrix
	
	
	for Inod=1:nod
		for Jnod=1:nod
			K=K+sparse(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
			
		end
	end
	
	%norm(max(kv-kv.'))
	
	K=(K+K.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	% Note: for numerical verificatin of distributed parameter gradient it is important to
	% not to use the complex conjugate transpose.
	
	
end



