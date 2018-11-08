function [beta2,Dbeta2Du] = calcBeta2in1d(u,C,m)
	
	% in the non-linear case beta2 goes to infty for u -> 0
	% a possibility of setting an upper limit to beta2 is by defining a minimum absolut velocity
    
	umin=min([1 mean(abs(u))/1000]);
	
	beta2= C.^(-1/m).* (abs(u)+umin).^(1/m-1) ;
	Dbeta2Du=(1/m-1)*C.^(-1/m).*(abs(u)+umin).^(1/m-2).*sign(u);

	
end

