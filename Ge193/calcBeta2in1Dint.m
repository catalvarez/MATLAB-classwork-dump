function [beta2int,Dbeta2Duuint] = calcBeta2in1Dint(uint,Cint,m)
    
    % calculated beta^2 and D beta^2/Du at integration point Iint
    
    SpeedZero=0.01;
    beta2int= Cint.^(-1/m).*(sqrt(uint.*uint+SpeedZero^2)).^(1/m-1) ;
    
    if nargout>1
        
        % The directional derivative is
        % D beta^2(u,v)[Delta u, Delta v]= (1/m-1) C^(-1/m) (u^2+v^2)^((1-3m)/2m)  (u \Delta u + v \Delta v)
        
        Dbeta2int=(1/m-1).*Cint.^(-1/m).*(uint.^2+SpeedZero^2).^((1-3*m)/(2*m));
        Dbeta2Duuint=Dbeta2int.*uint.*uint;
        
        
    end
end