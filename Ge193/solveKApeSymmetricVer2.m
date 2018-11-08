function  [x,y]=solveKApeSymmetricVer2(A,B,f,g)
	

	% x : solution vector
	% y : Lagrange parameters

	% Solves:
	%
	%  [A   B'] [x]= [f]
	%  [B  -C ] [y]  [g]
	%
	% where A is nxn, C is mxm , B is mxn
    % The system must be symmetric

    n=size(A,1) ; m=size(B,1);
    C=sparse(m,m);
    AA=[A B' ;B -C] ; bb=[f;g];
    sol=AA\bb;
    x=sol(1:n) ; y=sol(n+1:n+m);
    
    
end