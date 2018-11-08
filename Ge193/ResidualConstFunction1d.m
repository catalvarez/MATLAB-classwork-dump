function r=ResidualConstFunction1d(gamma,R,F,L,lambda,dlambda)
         
	
     if ~isempty(lambda)
        R=R+L'*(lambda+gamma*dlambda);
     end
    
    r=full(sqrt(R'*R/(F'*F)));
    
end

