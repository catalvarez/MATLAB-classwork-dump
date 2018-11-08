
function [Status,u,lambda,kv,rh,iteration]=...
        SSTREAM1dPicard(s,h,u,coordinates,connectivity,nip,etaInt,gfele,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,nPicard,tol,InfoLevel)
    
    Status=1;
    
    % Solves the SSTREAM equations in 1D (sparse and vectorized version)
    % derived from SSTREAM2d
    % lambda are the Lagrang parameters used to enfore the boundary conditions
    
    if any(h<0) ; error(' thickness negative ') ; end
    
    
    if n==1 && m==1 ;  % linear case
        %[beta2] = calcBeta2in1d(u,C,m);
        beta2= C.^(-1);
        [kv,rh]=kvrh1dPicard(s,h,coordinates,connectivity,nip,etaInt,gfele,beta2,alpha,rho,rhow,g) ;

        [u,lambda]=solveKApeSymmetricVer2(kv,L,rh,Lb);
        Status=0; iteration=0;
    else
        if InfoLevel > 1 ; disp(' non-linear ') ; end
        diffVector=zeros(100,1); diff=1e10;
        iteration_max=nPicard; iteration_min=2 ; nonlinear=1;
        iteration=0;
        while (diff > tol  && iteration < iteration_max )  || iteration < iteration_min
            iteration=iteration+1;
            
            [beta2] = calcBeta2in1d(u,C,m);
            %beta2= C.^(-1/m).* (sqrt(u.*u)).^(1/m-1) ;
            
            [kv,rh]=kvrh1dPicard(s,h,coordinates,connectivity,nip,etaInt,gfele,beta2,alpha,rho,rhow,g);
            
            
            ulast=u  ;
            [u,lambda]=solveKApeSymmetricVer2(kv,L,rh,Lb);

            D=mean(sqrt(u.*u));
            diffu=norm(u-ulast)/sqrt(length(u));
            diff=diffu/D ;               % normalized by the mean speed
            diffVector(iteration)=diff;
            
            if iteration>1 && InfoLevel> 1 ;disp([' Picard # ',num2str(iteration),' diff : ',num2str(diff)])  ;end
            [etaInt]=calcStrainRatesEtaInt1D(u,coordinates,connectivity,nip,AGlen,n);
        end
        iteration=iteration-1;    if diff <  tol ; Status=0; end
    end
    
    if InfoLevel>5 && nonlinear

            [detrended,a0,a1]=detrend_xt(log10(diffVector(1:iteration)),1:iteration);
            fprintf(' slope Picard  %f\n',a1)
            figure(700)
            semilogy(diffVector(1:iteration),'o-r') ; title('Picard') ;
            
            figure(800) ; hold off
            semilogy(diffVector(1:iteration),'o-r') ; title('Picard') ;
            
            hold on ; semilogy(10.^(a0+a1*[1:iteration]),'ko') ; hold off
    end
    
    
    
    
end