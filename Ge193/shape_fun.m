function [fun]=shape_fun(Iint,ndim,nod,points)  % 

% Iint   : integration point
% ndim   : dimensionality
% nip    : number of integration points
% nod    : number of nodes
% points : coordinates of integration points 

fun=zeros(nod,1) ;

switch ndim
  case 1
  %disp(' 1 d ')
  xi=points(Iint,1);

  switch nod
   case 2
      %disp(' 2 nodes ')
      fun(1)=(1-xi)/2 ;
      fun(2)=(1+xi)/2;
   case 3
     t1=-1-xi ;
     t2=-xi ;
     t3=1-xi;
     fun(1)=t2*t3/2 ;
     fun(2)=-t1*t3 ;
     fun(3)=t1*t2/2;
   case 4
     t1=-1-xi ;
     t2=-1/3-xi ;
     t3=1/3-xi ;
     t4=1-xi;
     fun(1)=t2*t3*t4*9/16  ;
     fun(2)=-t1*t3*t4*27/16;
     fun(3)=t1*t2*t4*27/16 ;
     fun(4)=-t1*t2*t3*9/16;
   case 5
     t1=-1 -xi ;
     t2=-0.5-xi ;
     t3=-xi ;
     t4=0.5-xi ;
     t5=1-xi;
     fun(1)=t2*t3*t4*t5*2/3 ;
     fun(2)=-t1*t3*t4*t5*8/3;
     fun(3)=t1*t2*t4*t5*d4 ;
     fun(4)=-t1*t2*t3*t5*8/3;
     fun(5)=t1*t2*t3*t4*2/3;
  otherwise
      disp(' not yet implemented ')
  end
otherwise
 disp(' not yet implemented ')
end