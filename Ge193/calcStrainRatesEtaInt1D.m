function [etaInt,xint,exx,DetaDuInt]=calcStrainRatesEtaInt1D(u,coordinates,connectivity,nip,AGlen,n)
    
    
    %%
    
    %load TestSave u coordinates connectivity nip AGlen n
    
    % calculates strain rates and effective viscosity at integration points
    % vectorized
    %
    
    InfoLevel=0;
    
    e00=1e-5;  % min strain to create a `Null' viscosity
    
    
    [Nele,nod]=size(connectivity); ndim=1;
    %    [points,weights]=sample('triangle',nip,ndim);
    [points,weights]=sample(ndim,nip,nod);
    
    unod=reshape(u(connectivity,1),Nele,nod);
    
    
    coox=reshape(coordinates(connectivity,1),Nele,nod);
    
    exx=zeros(Nele,nip);
    xint=zeros(Nele,nip) ;
    
    for Iint=1:nip                           % loop over integration points
        fun=shape_fun(Iint,ndim,nod,points) ; % nod x 1   : [N1 ; N2 ; N3] values of form functions at integration points
        
        
        xint(:,Iint)=coox*fun;
        
        [Deriv]=derivVector1d(coordinates,connectivity,nip,Iint); % Nele x dof x nod
        
        for I=1:nod
            exx(:,Iint)=exx(:,Iint)+Deriv(:,1,I).*unod(:,I);
        end
        
    end
      
    e=abs(exx)+e00;
    
    %fprintf(' max(abs(exx(:))) min(abs(exx(:))) %g %g \n ',max(abs(exx(:))),min(abs(exx(:))))
    
    etaInt=0.5*AGlen.^(-1/n).*e.^((1-n)/n);
    DetaDuInt=(1-n).*AGlen.^(-1/n).*e.^(1/n-2).*sign(exx)/ 2/ n ;
      
    if InfoLevel>1
        figure ; plot(xint(:),etaInt(:),'ro') ; title('etaInt')
        figure ; plot(xint(:),DetaDuInt(:),'ro') ; title('DetaDuInt')
        figure ; plot(xint(:),exx(:),'ro') ; title('exx')
        input(' Press return to proceed ')
    end
    
    %%
    
end
