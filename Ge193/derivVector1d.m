function [Deriv,detJ]=derivVector1d(coordinates,connectivity,nip,Iint)
    
    % calculates the derivatives of form functions with respect to x and y and the Jacobian
    % at a given integration point (This function must be called within a loop over integration points)
    
    % Deriv : Nele x dof x nod
    %  detJ : Nele
   
    [Nele,nod]=size(connectivity) ; ndim=1; dof=1;
    
	
    %[points,weights]=sample('triangle',nip,ndim);
	[points,weights]=sample(ndim,nip,nod);
	
    Deriv=zeros(Nele,dof,nod);
    
    %hnod=reshape(h(connectivity,1),Nele,nod);
    coox=reshape(coordinates(connectivity,1),Nele,nod);
         
   % fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
    der=shape_der(Iint,ndim,nod,points);  % dof x nod : dNj/dXi=[dN1/dx dN2/dx dN3/dx; dN1/dy dN2/dy dN3/dy]
    %hint=hnod*fun;
    
    % calculate J, detJ, and deriv in vectorized way
    %  x=x_p N_p   ; h=h_q N_q
    % der=dN_p/dxi_q
    % J11=dNp/dxi_1 x_p
    % (x,y) and N(r,s)
    % x=x_p N_p(r,s) , y=y_p N_p(r,s) , h=h_p N_p(r,s)
    %
    % dh/dx=dh/dr dr/dx+ dh/ds ds/dx =  h_p (dN_p/dr dr/dx + dN_p/ds ds/dx)
    %
    %  der(1:dof,1:nod)=[ dN1/dr  dN2/dr ... dNnod/dr ]
    %                   [ dN1/ds  dN2/ds ... dNnod/ds ]
    %
    % The Jakobian is: J(1:dof,1:dof)= [dx/dr  dy/dr ]
    %                                  [dx/ds  dy/ds ]
    % and can be calculated as
    % J=der*coo; % (dof x nod) x (nod x dof) = dof x dof
    % this expression is only evaluated at the integration points,
    % der is the derivative of from functions evaluated at integration points
    %
    % At each integration point
    % J11(1:Nele)=dx/dr=x_p(1:Nele) dN_p/dr = coox*der(1,:)'
    % where coox(1:Nele,1:nod)
    % J12(1:Nele)=dy/dr=cooy*der(1,:)'
    % deriv=inv(J)*der=iJ(1:dof,1:dof) der(1:dof,1:nod) = dof x nod
    %
    % deriv=[ dN1/dx dN2/dx ... dNnod/dx ] = [dr/dx ds/dx] [dN1/dr   dN2/dr ... dNnod/dr ]
    %       [ dN1/dy dN2/dy ... dNnod/dy ]   [dr/dy ds/dy] [dN1/ds   dN2/ds ... dNnod/ds ]
    %
    % iJac=[dr/dx ds/dx]
    %      [dr/dy ds/dy]
    %
    % iJac=[J22 -J12]  / detJ
    %      [-J21 J22]
    %
    % dN1/dx=iJ11*der(1,1)+iJ12*der(2,1)
    % dN2/dx=iJ11*der(1,2)+iJ12*der(2,2)
    % dNnod/dx=iJ11*der(1,nod)+iJ12*der(2,nod)
    % dNnod/dy=iJ21*der(2,nod)+iJ22*der(2,nod)
    % deriv=iJ der; % (dof x dof) x (dof x nod) = dof x nod
    %
     
    % this is the Jacobian 
    J11=coox*der(1,:)' ; % dN1/dxi1 : (Nele x nod) x nod = Nele
    
    detJ=J11 ; % This is the determinant of the Jacobian at this integration point
       
    for J=1:nod
    
	    Deriv(:,1,J)=der(1,J)./J11;
    end
    
    
end




