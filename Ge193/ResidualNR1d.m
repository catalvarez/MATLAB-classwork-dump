function r=ResidualNR1d(gamma,F,L,s,h,u,lambda,du,dlambda,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g)
    
    R=rh1NR(s,h,u+gamma*du,AGlen,n,C,m,coordinates,connectivity,nip,gfele,alpha,rho,rhow,g);

    r=ResidualConstFunction1d(gamma,R,F,L,lambda,dlambda);
    
   
    
    
end