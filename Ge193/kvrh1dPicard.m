
function [kv,rh]=kvrh1dPicard(s,h,coordinates,connectivity,nip,etaInt,gfele,beta2,alpha,rho,rhow,g)
	
	% calculates the FE matrix and right-hand side in a vectorized form
	
	Nnodes=max(connectivity(:)); [Nele,nod]=size(connectivity); 
	ndim=1; dof=1; neq=dof*Nnodes;
	
	
	hnod=reshape(h(connectivity,1),Nele,nod);   % Nele x nod
	snod=reshape(s(connectivity,1),Nele,nod);
	
	
	
	

	beta2nod=reshape(beta2(connectivity,1),Nele,nod);
	
	
	rhog=rho*g; ca=cos(alpha); sa=sin(alpha);
	
	
	%[points,weights]=sample('triangle',nip,ndim);
	[points,weights]=sample(ndim,nip,nod);
	
	
	% get local coordinates and weights
	
	
	kv=sparse(neq,neq);
	d1d1=zeros(Nele,nod,nod);
	b1=zeros(Nele,nod);
	
	for Iint=1:nip
		
		
		fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
		[Deriv,detJ]=derivVector1d(coordinates,connectivity,nip,Iint);
		
		
		
		% Deriv : Nele x dof x nod
		%  detJ : Nele
		
		% values at integration this point
        
        % gfint=gfnod*fun;
        % gfint=gfele
        
        
		hint=hnod*fun;
		%sint=snod*fun;
		%uint=unod*fun;
		%vint=vnod*fun;
		beta2int=beta2nod*fun; 
		 
		
		%beta2int=beta2int.*gfint(:,Iint);
       
        beta2int=beta2int.*gfele;  % all or none integration points of a element grounded
        
		etaint=etaInt(:,Iint) ;  % I could consider calculating this here
		
		 		
		
		% derivatives at this integration point for all elements
		dsdx=zeros(Nele,1); dhdx=zeros(Nele,1);
		
		for Inod=1:nod
			
			dsdx=dsdx+Deriv(:,1,Inod).*snod(:,Inod);
			dhdx=dhdx+Deriv(:,1,Inod).*hnod(:,Inod);
			
			
		end
		
				
		
		detJw=detJ*weights(Iint);
		
		
		for Inod=1:nod
			for Jnod=1:nod
				
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+4*etaint.*hint.*Deriv(:,1,Inod).*Deriv(:,1,Jnod).*detJw;
				d1d1(:,Inod,Jnod)=d1d1(:,Inod,Jnod)+beta2int.*fun(Jnod).*fun(Inod).*detJw;
								
			end
			
			t1=rhog*hint.*((dsdx-(1-rho/rhow).*dhdx)*ca-sa).*fun(Inod);
			t2=-0.5*ca*rhog.*(1-rho/rhow).*hint.^2.*Deriv(:,1,Inod); % correction on 11 march 2010
		    
			b1(:,Inod)=b1(:,Inod)-(t1+t2).*detJw;
			
			
			
		end
	end
	
	
	
		
	% assemble right-hand side
	
	rh=sparse(neq,1);
	for Inod=1:nod
		
		rh=rh+sparse(connectivity(:,Inod),ones(Nele,1),b1(:,Inod),neq,1);
		
		
	end
	
	rh=full(rh);
	
	% assemble matrix
	
	
	for Inod=1:nod
		for Jnod=1:nod
			kv=kv+sparse(connectivity(:,Inod),connectivity(:,Jnod),d1d1(:,Inod,Jnod),neq,neq);
		
		end
	end
	
	
	
	kv=(kv+kv.')/2 ; % I know that the matrix must be symmetric, but numerically this may not be strickly so
	% Note: for numerical verificatin of distributed parameter gradient it is important to
	% not to use the complex conjugate transpose.
	
	
end



