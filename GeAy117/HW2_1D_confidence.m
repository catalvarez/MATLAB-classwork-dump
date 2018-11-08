function [ lowerbound,upperbound,peak ] = HW2_1D_confidence( x,Px,interval )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
dx = x(2)-x(1);

Px = Px./sum(Px.*dx);

[Ppeak, index] = max(Px);
peak = x(index);
cdf = zeros(length(Px),1);

cdf(1) = Px(1)*dx;
for i = drange(2:length(Px))
    cdf(i) = cdf(i-1)+Px(i)*dx;
end

    ub = 0.5+0.005*interval;
    lb = 0.5-0.005*interval;

    lindex = find(cdf<ub,1,'last');
    uindex = lindex+1;
upperbound = (x(uindex)-x(lindex))/(cdf(uindex)-cdf(lindex))*(ub-cdf(uindex))+x(uindex);
    lindex = find(cdf<lb,1,'last');
    uindex = lindex+1;
lowerbound = (x(uindex)-x(lindex))/(cdf(uindex)-cdf(lindex))*(lb-cdf(uindex))+x(uindex);

end

