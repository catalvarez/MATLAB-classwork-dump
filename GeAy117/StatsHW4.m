function f = StatsHW4(a)

% t = load('t.mat');
% Y = load('Y.mat');
% sigma = load('sigma.mat');
mat = dlmread('ps4_data.txt');
f = 0.5.*sum((a(1).*sin(2*pi*mat(:,1)./a(3))+a(2)-mat(:,2)).^2./mat(:,3).^2);

end