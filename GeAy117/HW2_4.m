alphareal = 6;
betareal = 3;

% thetapos = rand(128,1).*pi-pi./2;
% lighthouse2 = betareal.*tan(thetapos)+alphareal;

da = 0.1;
alpha = -10:da:10;
db = 0.1;
beta = 0:db:10;
L = zeros(length(alpha),length(beta));

for i = drange(1:length(alpha))
    for j = drange(1:length(beta))
L(i,j) = sum(log(beta(j)./pi./(beta(j).^2+(lighthouse2-alpha(i)).^2)));
    end
end

pdf2D = L-max(max(L))+1;

P = exp(pdf2D)./sum(sum(exp(pdf2D)*da*db));
figure (1)
contour(beta,alpha,P)
% xlim([0 10])
% ylim([-10 10])
title('The Lighthouse')
xlabel('Distance from Shore (km)')
ylabel('Location on Shoreline (km)')

B = sum(P*da,1);
B = B/sum(B*db);
A = sum(P*db,2);
A = A/sum(A*da);
figure (2)
plot(beta,B)
title('The Lighthouse: Beta')
xlabel('Distance from Shore (km)')
ylabel('Probability')
figure (3)
plot(alpha,A)
title('The Lighthouse: Alpha')
xlabel('Location on Shoreline (km)')
ylabel('Probability')

[ alowerbound,aupperbound,apeak ] = HW2_1D_confidence(alpha,A,68);
[ blowerbound,bupperbound,bpeak ] = HW2_1D_confidence(beta,B,68);
