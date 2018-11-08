function [u,lambda,kv,rh]=SSTREAM1d(s,h,u,coordinates,connectivity,nip,gfele,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,NumConstants)

    Status=1;
    tryNR=1; 

    if n==1 && m==1 ; tryNR=0 ; end
    
    itNR=0 ; itPicard=0; IT=0;
    
    while Status==1 && IT<50
       
        if tryNR==1
            [Status,u,lambda,kv,R,iterationNR]=...
                SSTREAM1dNR(s,h,u,coordinates,connectivity,nip,gfele,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,NumConstants);
                rh=-R;
                
        else
            Status=1; iterationNR=0;
        end
        
        itNR=itNR+iterationNR;
        
        if Status ~= 0
            if NumConstants.InfoLevel > 1 ; disp(' NR did not converge, switch to Picard iteratoin ') ; end
            
            [etaInt]=calcStrainRatesEtaInt1D(u,coordinates,connectivity,nip,AGlen,n);
            
            if tryNR ; end
            [Status,u,lambda,kv,rh,iterationPicard]=...
                SSTREAM1dPicard(s,h,u,coordinates,connectivity,nip,etaInt,gfele,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,NumConstants.MaxPicard,NumConstants.Picardtol,NumConstants.InfoLevel);
            
            itPicard=itPicard+iterationPicard;
        end
        
        IT=itPicard+itNR;
    end
    
    
    if NumConstants.InfoLevel>0
        fprintf('iterationNR %i iterationPicard %i \n',itNR, itPicard)
    end

    if Status==1 ; fprintf('Warning : non-linear iteration in SSTREAM1d did not fully converge \n') ; end
    
end
