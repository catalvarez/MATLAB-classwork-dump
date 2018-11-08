t = y+(m-1)./12;
for i = 1:length(d)
    a = strcmp('nan',d(i));
    if a == 1
        d(i) = 0
    end
end
        
plot(t,d)