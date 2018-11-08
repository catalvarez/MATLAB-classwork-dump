R = 20000;
L = 30000;
x = -L:100:L;
y = -L:100:L;
q = zeros(2*L/100,2*L/100);
B = 10^-11;

% for i = 1:2*L/100
%     for j = 1:2*L/100
%         r = sqrt(x(j)^2+y(i)^2);
%         I = r < R;
%         if I == 1
%           q(i,j) = -B/(4*R).*((y(i)-2*R).^2.-5*R^2.+x(j)^2);
%         else
%           q(i,j) =  B*y(i);
%         end
%     end
% end

for i = 1:2*L/100
    for j = 1:2*L/100
        r = sqrt(x(j)^2+y(i)^2);
        I = r < R;
        if I == 1
          q(i,j) = -2*B/R.*((y(i)-.25*R).^2.-(.25*R)^2.+x(j)^2-R^2);
        else
          q(i,j) =  B*y(i);
        end
    end
end



contour(q)
title('y0 = 0.25R')
xlabel('Longitude')
ylabel('Latitude')
    