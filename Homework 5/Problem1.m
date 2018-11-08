phos = [0.1:0.1:10];
sil = [1:100];

%phosphate limiting (silicate replete)
% y1 = 0.27.*phos./(0.03+phos);
% plot(phos,y1)
% hold
% y2 = 0.25.*phos./(0.03+phos);
% h = plot(phos,y2)
% set(h, 'Color', 'Green')
% xlabel('Concentration')
% ylabel('Rate')
% title('Phosphate Limiting')

%silicate limiting
% y1 = 0.27.*sil./(7+sil);
% plot(sil,y1)
% hold
% y2 = 0.25.*sil./(sil);
% h = plot(sil,y2)
% set(h, 'Color', 'Green')
% xlabel('Concentration')
% ylabel('Rate')
% title('Silicate Limiting')

sil2 = sil';
diatom1 = 0.27.*(sil2./(7+sil2))*(phos./(phos));
diatom2 = 0.27.*(sil2./sil2)*(phos./(0.03+phos));
diatom = zeros(100,100);
for i = [1:100];
    for j = [1:100];
        diatom(i,j) = min(diatom1(i,j),diatom2(i,j));
    end
end
    
other = 0.25.*(sil2./sil2)*(phos./(0.03+phos));
ratio = diatom./other;

% [C,h] = contourf(phos,sil2,diatom)
%  hold
% xlabel('Phosphate Concentration')
% ylabel('Silicate Concentration')
% title('Diatom Growth Rate')

% [C,h] = contourf(phos,sil2,other)
% xlabel('Phosphate Concentration')
% ylabel('Silicate Concentration')
% title('Phytoplankton Growth Rate')

[C,h] = contourf(phos,sil2,ratio)
clabel(C,h)
xlabel('Phosphate Concentration')
ylabel('Silicate Concentration')
title('Diatom/Phytoplankton Growth Rate')
