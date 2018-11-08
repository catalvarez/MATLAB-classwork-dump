function [points,weights]=sample(ndim,nip,nod);  % get local coordinates and weights

% ndim : dimensionality
% nip  : number of integration points
% nod  : number of nodes

% output: coordinates of integration points and corresponding weights for numerical integration

points=zeros(nip,ndim) ; weights=zeros(nip,1);

switch ndim
  case 1
  %disp(' 1 d ')
  switch nip
    case 1
     points(1,1)=0 ;
     weights(1)=2 ;
    case 2
     %disp(' two integration points ')
     points(1,1)=-0.577350269189626;
     points(2,1)= 0.577350269189626;
     weights(1) = 1.000000000000000;
     weights(2) = 1.000000000000000;
    case 3
     %disp(' three integration points ')
     points(1,1)=-0.774596669241484;
     points(2,1)= 0.000000000000000;
     points(3,1)= 0.774596669241484;
     weights(1) = 0.555555555555556;
     weights(2) = 0.888888888888889;
     weights(3) = 0.555555555555556;
   case 4 
     points(1,1)=-0.861136311594053;
     points(2,1)=-0.339981043584856;
     points(3,1)= 0.339981043584856;
     points(4,1)= 0.861136311594053;
     weights(1) = 0.347854845137454;
     weights(2) = 0.652145154862546;
     weights(3) = 0.652145154862546;
     weights(4) = 0.347854845137454;
   case 5
     points(1,1)=-0.906179845938664;
     points(2,1)=-0.538469310105683;
     points(3,1)= 0.000000000000000;
     points(4,1)= 0.538469310105683;
     points(5,1)= 0.906179845938664;
     weights(1) = 0.236926885056189;
     weights(2) = 0.478628670499366;
     weights(3) = 0.568888888888889;
     weights(4) = 0.478628670499366;
     weights(5) = 0.236926885056189;
   case 6
     points(1,1)=-0.932469514203152;
     points(2,1)=-0.661209386466265;
     points(3,1)=-0.238619186083197;
     points(4,1)= 0.238619186083197;
     points(5,1)= 0.661209386466265;
     points(6,1)= 0.932469514203152;
     weights(1) = 0.171324492379170;
     weights(2) = 0.360761573048139;
     weights(3) = 0.467913934572691;
     weights(4) = 0.467913934572691;
     weights(5) = 0.360761573048139;
     weights(6) = 0.171324492379170;
  otherwise
      disp(' not yet implemented ')
  end
otherwise
 disp(' not yet implemented ')
end


 