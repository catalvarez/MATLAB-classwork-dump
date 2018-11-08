function [der]=shape_der(Iint,ndim,nod,points)  % 

% Iint : integration point
% ndim : dimensionality
% nip  : number of integration points
% nod  : number of nodes
   zero=0.0; pt125=0.125; pt25=0.25; pt5=0.5;   
   pt75=0.75; one=1.0; two=2.0; d3=3.0; d4=4.0; d5=5.0; 
   d6=6.0; d8=8.0; d9=9.0; d10=10.0; d11=11.0;             
   d12=12.0; d16=16.0; d18=18.0; d27=27.0; d32=32.0;       
   d36=36.0; d54=54.0; d64=64.0; d128=128.0;

switch ndim
  case 1
  %disp(' 1 d ')
  xi=points(Iint,1);
  switch nod
      case 2
      %disp(' two nodes ')
      der(1,1)=-0.5 ;
      der(1,2)=0.5 ;
   case 3
     t1=-one-xi ;
     t2=-xi  ;
     t3=one-xi;
     der(1,1)=-(t3+t2)/two  ;
     der(1,2)=(t3+t1)    ;
     der(1,3)=-(t2+t1)/two   ;
   case 4
     t1=-one-xi ;
     t2=-one/d3-xi ;
     t3=one/d3-xi ;
     t4=one-xi;
     der(1,1)=-(t3*t4+t2*t4+t2*t3)*d9/d16     ;
     der(1,2)=(t3*t4+t1*t4+t1*t3)*d27/d16 ;
     der(1,3)=-(t2*t4+t1*t4+t1*t2)*d27/d16 ;
     der(1,4)=(t2*t3+t1*t3+t1*t2)*d9/d16   ;
   case 5
     t1=-one-xi ;
     t2=-pt5-xi ;
     t3=-xi ;
     t4=pt5-xi ;
     t5=one-xi;
     der(1,1)=-(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3   ;
     der(1,2)=(t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3;
     der(1,3)=-(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4 ;
     der(1,4)=(t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3;
     der(1,5)=-(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3;
  otherwise
      disp(' not yet implemented ');
  end;
otherwise;
 disp(' not yet implemented ');
end;
