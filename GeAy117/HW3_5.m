da = 0.001;
alpha = 0:da:.1;
db = 0.1;
beta = .5:db:2.5;
dm = 0.01;
m = -2:dm:2;
L = zeros(length(alpha),length(beta));
for k = drange(1:length(m))
    trial = exp(-(m(k)-M).^2./(2*0.04).^2);
end

Mi = sum(trial);

for i = drange(1:length(alpha))
    for j = drange(1:length(beta)) 
%         for k = drange(1:length(m))
L(i,j) = prod((alpha(i).*10.^(beta(j)*Mi)).^d.*(1-(alpha(i).*10.^(beta(j)*Mi))).^(1-d));%.*exp(-(m(k)-M).^2./(2*0.04).^2)));
%         end
    end
end

% L = exp(L);
P = L./sum(sum(L*da*db));
figure (1)
contour(beta,alpha,P)
colorbar
% xlim([0 10])
% ylim([-10 10])
title('Stars with Binary Companions: Metallicity')
xlabel('Beta value')
ylabel('Alpha value')

B = sum(P*da,1);
B = B/sum(B*db);
A = sum(P*db,2);
A = A/sum(A*da);
figure (2)
plot(beta,B)
title('Stars: Beta')
xlabel('Beta Value')
ylabel('Probability')
figure (3)
plot(alpha,A)
title('Stars: Alpha')
xlabel('Alpha Value')
ylabel('Probability')

[ alowerbound,aupperbound,apeak ] = HW2_1D_confidence(alpha,A,68);
[ blowerbound,bupperbound,bpeak ] = HW2_1D_confidence(beta,B,68);
