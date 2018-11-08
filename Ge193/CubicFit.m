function [xmin] = CubicFit(Slope,g0,g1,g2,l1,l2)
    
    % finds the minimum based on a cubic fit given three function values
    % a 0, l1, and l2 and the slope at 0, i.e.
    % g0=f(0)  ; g1=f(l1) ; g2=f(l2) with Slope=f'(0) 
    
    ab=([1/l1^2 -1/l2^2 ; -l2/l1^2 l1/l2^2 ]* [g1-Slope*l1-g0 ; g2-Slope*l2-g0])/(l1-l2);
    xmin=(-ab(2)+sqrt(ab(2)^2-3*ab(1)*Slope))/3/ab(1);
    
    if ab(1) == 0 
        figure; plot([0,l1,l2],[g0,g1,g2])
        hold on
        plot([0,l1*0.01],[g0,g0+0.01*Slope],'g')
        Slope
        error('sdaf')
    end
    
end

