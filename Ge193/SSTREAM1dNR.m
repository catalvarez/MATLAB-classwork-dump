
function [Status,u,lambda,K,R,iteration]=...
        SSTREAM1dNR(s,h,u,coordinates,connectivity,nip,gfele,AGlen,C,L,Lb,lambda,n,m,alpha,rho,rhow,g,NumConstants)
 
    % Solves the SSTREAM equations in 1D (sparse and vectorized version)
    % Uses Newton-Raphson to solve the non-linear equations
    % derived from SSTREAM2d
    % lambda are the Lagrange parameters used to enfore the boundary conditions
    
    % I'm solving kv(u_k) du_k =rh(r_k)
    %  with   u_{k+1}= u_k+ du_k
    % the solutino has converged once du < tol
    % another measure of convergence is the size of rh which must go to zero
    % I use du/mean(d) as a measure for convergence
    
    %Line search
    % write u=u+ gamma du and require
    %  0=< du , rh(u+gamma du) >  = du'*rh(u+gamma du)  =:R(gamma)
    % where kv du = rh(u)
    %
    %
    % Note that the derivative at gamma=0 is : dR/dgamma=-du' kv du = -du rh(u)=R(0)
    % if I know the values at R(0) and R(1) and the derivative at gamma=0 is R(0) it follows that
    % R(gamma)=(1-gamma) R(0) + R(1) gamma^2 = R(0)- gamma R(0) + R(1) gamma^2
    %
    %  The minumun is found for : gamma= (R(0) pm sqrt( R(0)^2-4 R(0) R(1)))/(2 R(1))
    
    %
    % R(0)=du'*rh(u)  and R(1)=du'*rh(u+du);
    
     
    tol=NumConstants.NRtol;
    InfoLevel=NumConstants.InfoLevel;
    nNR=NumConstants.MaxNR;
    
    if any(h<0) ; error(' thickness negative ') ; end
    
      Status=0 ;  % Returns 0 if everything seems OK and it converged to prescribed tolerance
    
     dlambda=zeros(length(lambda),1);
    
    
    
    
    diffVector=zeros(100,1);  diffr=1e10;
    
    
    
    if n==1 && m==1 ;
        iteration_max=1 ; iteration_min=1 ;
    else
        if InfoLevel > 2 ; disp(' non-linear ') ; end
        iteration_max=nNR; iteration_min=2 ;
    end
    
    
    
    iteration=0;
    while (diffr > tol  && iteration < iteration_max )  || iteration < iteration_min
        iteration=iteration+1;
        
        
        [K,R,T,F0]=kvrh1NR(s,h,u,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);
                
        r0=ResidualConstFunction1d(0,R,F0,L,lambda,dlambda);
        
        [du,dlambda]=solveKApeSymmetricVer2(K,L,-R-L'*lambda,Lb-L*u);
        Slope=-r0;
        
        
        R=rh1NR(s,h,u+du,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);
        r1=ResidualConstFunction1d(1,R,F0,L,lambda,dlambda);
        
        
        
        gamma=1 ; rgamma=r1;
        
        if rgamma/r0 > 1  % overshooting!
            
            if InfoLevel>3
                fprintf(' Newton step overshot with r1/r0 %g so try backtracking \n',r1/r0)
            end
            
            gamma=-Slope/2/(r1-r0-Slope);  % location of minimum based on a quadradic fit
            gamma=max([0.1 min([0.8 gamma])]);
            
            % just go directly to a cubic fit without trying using the value obtained by a quadric fit
            
            b=gamma; rb=ResidualNR1d(b,F0,L,s,h,u,lambda,du,dlambda,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);

            
            % and find minimum based on a cubic fit
            [gamma] = CubicFit(Slope,r0,rb,r1,b,1);
            
            gamma=max([0.1*b min([0.5*b gamma])]);
            
            rgamma=ResidualNR1d(gamma,F0,L,s,h,u,lambda,du,dlambda,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);
            
            
            rvalues=[r1 ; rb ; rgamma] ; gammas=[1 ; b ; gamma];
            [rgamma,I]=min(rvalues) ; gamma=gammas(I);
            if InfoLevel>3
                fprintf(' After one cubic backtracking gamma is %g with a ratio %g \n',gamma,rgamma/r0)
            end
        end
        
        if rgamma/r0 > 0.8  % try line search
            if InfoLevel>3
                fprintf(' try line search as ratio is %g \n ',rgamma/r0)
            end
            gammaOld=gamma ; rgammaOld=rgamma;
           
            
            [gamma,rgamma]=fminbndGHG(...
                @(gamma) ResidualNR1d(gamma,F0,L,s,h,u,lambda,du,dlambda,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g),...
                0.01,10,gamma,rgamma,optimset('TolX',0.01,'Display','off','MaxFunEvals',10));
            if InfoLevel>3
                fprintf('Line search gives gamma %g and rgamma %g with the ratio %g \n',gamma,rgamma,rgamma/r0)
            end
            
            if rgamma> rgammaOld ; rgamma=rgammaOld ; gamma=gammaOld ; end
            
            if InfoLevel>10
                NN=40 ; gammavector=zeros(NN,1) ; rvector=zeros(NN,1);
                for I=1:NN
                    gammavector(I)=2*I/NN;
                    rvector(I)=ResidualNR1d(gammavector(I),F0,L,s,h,u,lambda,du,dlambda,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);
                end
                
                figure(999), plot(gammavector,rvector) ;
                hold on ; plot(gamma,rgamma,'o') ; input('Return when ready')  ;  hold off
            end
            
        end
        
        if rgamma/r0>0.999 ;  Status=1 ; break ; end
        
        
        u=u+gamma*du;
        lambda=lambda+gamma*dlambda;
        
        D=mean(sqrt(u.*u));
        diffu=norm(du)/sqrt(length(du));
        diff=diffu/D ;               % normalized by the mean speed
        
        
        diffr=rgamma ;
        diffVector(iteration)=rgamma;
        
        if InfoLevel> 1 ; fprintf(' NR # %g Residual %g uDiff %g gamma %g ratio %g \n',iteration,rgamma,diff,gamma,rgamma/r0) ; end
        
        
    end
    
    if InfoLevel>10 && iteration >= 2
        
        
        % 		%lndu=c0+c1 * iteration
        % 		N=3;
        % 		[detrended,a0,a1]=detrend_xt(log10(diffVector(1:N)),1:N);
        % 		a1
        
        
        figure(2000) ;
        semilogy(diffVector(1:iteration),'x-g') ; title('NR') ;
        
        
        N=max([1,iteration-5]);
        
        [detrended,a0,a1]=detrend_xt(log10(diffVector(N:iteration)),N:iteration);
        fprintf(' slope NR : %g \n',a1)
        
        figure(2100) ;
        semilogy(diffVector(1:iteration),'x-g') ; title('NR') ;
        hold on ; semilogy(10.^(a0+a1*[1:iteration]),'ko') ;
        
    end
    
    % Check if iteration converged within max number of steps allowed
    % but also don't consider it a failure if residual < 1e-10  
    if iteration == iteration_max && diffr > 1e-10 
        Status=1;
        if InfoLevel>0 ; disp(' maximum number of iteration reached in NR! ') ; end
        
    end
    
    iteration=iteration-1;
    
end


