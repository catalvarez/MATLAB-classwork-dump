%taking rows out using flag columns (2=good data)
% dataf1 = a09hy1(a09hy1(:,7)==2,:); 
% dataf2 = dataf1(dataf1(:,9)==2,:);
% dataf3 = dataf2(dataf2(:,11)==2,:);
% data = dataf3(dataf3(:,13)==2,:);
% pdepth = sw_dpth(data(:,4),data(:,1));
 
%problem 2a
% plot(data(:,6),data(:,5),'.')
% ylabel('Temperature [K]')
% xlabel('Salinity')

%problem 2b
% minlong = min(data(:,2));
% maxlong = max(data(:,2));
% maxdepth = -min(pdepth);
% mindepth = -max(pdepth);
% [xq,yq] = meshgrid(minlong:maxlong, mindepth:maxdepth);
% 
%change data(:,x) and title to plot each variable
% zonal = griddata(data(:,2),-pdepth,data(:,12),xq,yq);
% contourf(xq,yq,zonal,100,'linestyle','none')
% colorbar()
% title('Nitrate')
% xlabel('Longitude')
% ylabel('Depth')

%problem 2c
% pcadata = data;
% pcadata(:,13) = [];
% pcadata(:,11) = [];
% pcadata(:,9) = [];
% pcadata(:,7) = [];
% pcadata(:,4) = [];
% pcadata(:,3) = [];
% pcadata(:,2) = [];
% pcadata(:,1) = [];

% R = corrcoef(pcadata);
% [V,Lambda] = eigsort(R);
% AR = V*(Lambda.^0.5);
% PoV = 100*diag(Lambda)/trace(Lambda);
% disp([PoV cumsum(PoV)]);
% 
%problem 2d
% plot(AR)
% ylabel('Loading')
% set(gca,'xtick',[1 2 3 4 5]);
% set(gca,'xticklabel',{'Temperature';'Salinity';'Oxygen';'Silicate';'Nitrate'});