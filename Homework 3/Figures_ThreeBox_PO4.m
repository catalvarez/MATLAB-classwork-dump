concmatrix = [0.0000 0.0002 0.0026; 0.0000 0.0004 0.0026; 0.0000 0.0006 0.0026; 0.0000 0.0007 0.0026; 0.0000 0.0009 0.0026; 0.0000 0.0010 0.0026; 0.0000 0.0011 0.0026; 0.0000 0.0012 0.0026; 0.0000 0.0012 0.0026; 0.0000 0.0013 0.0026];
eff = 1.-concmatrix(:,2)./concmatrix(:,3);
hdflux = [4;8;12;16;20;24;28;32;36;40];
plot(hdflux,eff)
%hold
title('Productivity, First Order Rate Constant')
xlabel('f_H_D [Sv]')
ylabel('Efficiency')

concmatrix2 = [0.0000 0.0010 0.0075; 0.0000 0.0018 0.0075; 0.0000 0.0024 0.0075; 0.0000 0.0029 0.0075; 0.0000 0.0033 0.0075; 0.0000 0.0036 0.0075; 0.0000 0.0039 0.0075; 0.0000 0.0041 0.0075; 0.0000 0.0044 0.0075; 0.0000 0.0045 0.0075];
eff2 = 1.-concmatrix2(:,2)./concmatrix2(:,3);
hdflux = [4;8;12;16;20;24;28;32;36;40];
plot(hdflux,eff2)
title('Productivity, Constant Flux')
xlabel('f_H_D [Sv]')
ylabel('Efficiency')