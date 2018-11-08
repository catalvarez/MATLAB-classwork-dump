da = 0.1;
alpha = -10:da:10;
alphaguess = -1;
beta = 2;
L = zeros(length(alpha),1);
prior = ones(length(alpha),1);

for i = drange(1:length(alpha))
L(i) = sum(log(beta./pi./(beta.^2+(lighthouse1-alpha(i)).^2)));
prior(i) = 1./(2*pi).*exp(-(alpha(i)-alphaguess).^2./2);
end

L = L-min(L)+1;
figure (1)
plot(alpha,L)
pdf1D = L.*prior;

gaussian_P = exp(pdf1D)./sum(exp(pdf1D)*da);
uniformprior_P = exp(L)./sum(exp(L)*da);
figure(3)
plot(alpha,gaussian_P)
hold on
plot (alpha,uniformprior_P,'r')
title('The Lighthouse')
xlabel('Location on Shoreline (km)')
ylabel('Probability of Location')

[ lowerbound,upperbound,peak ] = HW2_1D_confidence(alpha,Prob,68);