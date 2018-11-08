function [trans]=SIS_Tsc_t_3d_m(kx,ky,t,alpha,H,eta,C,rho,g,m)

% SSTREAM

% time-dependent ratio between surface topography and
% basal slipperiness Fourier space in dimensional units

% for non-dimensional a la gudmundsson 2003
% put eta=1/2, H0=1,rho*g=1/sin(alpha), lambda in units of mean ice thickness
% Cnondimensional=2 eta/H0 Cdimesional, and then
% U0nondimensional=Cnondimensional 

if isvector(kx) | isvector(ky)
  % if k & l are vectors then 
  kx=kx(:); ky=ky(:);
  k=repmat(kx,1,length(ky));
  l=repmat(ky',length(kx),1);
else
  k=kx ; l=ky;
end


j2=k.^2+l.^2;
ca=cot(alpha);
tau=rho*g*sin(alpha)*H;
U=C*tau^m;
gamm=tau^(1-m) / (C*m); % can't use gamma because it is a function

p=i*k*U - ((k.*(-i + ca*H*k) + ca*H.*l.^2)*tau)./(gamm + 4*H*j2*eta);
xi=gamm+4*H*j2*eta;

if isnan(t)
  expp=0;
  disp(' steady state option used ')
else
  expp=exp(p*t);
end

t1=i*(expp-1)*H.*k*U*gamm;
t2=p.*xi;

trans=t1./t2;


return

end



