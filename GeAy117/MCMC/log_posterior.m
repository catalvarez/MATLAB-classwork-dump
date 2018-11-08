function lnL = log_posterior(param,X,Y,sigma)
% takes parameters in the form param(value,param#) and data and solves for the natural log of the
% likelihood of the data given the model

% Modify here for different model structures
% y = param(1).*X+param(2);
y = param(1).*cos(2.*pi./param(2).*X+param(3))+param(4);

    % solves lnL for normally distributed uncertainties
    lnL = -0.5.*sum((y-Y).^2./(sigma.^2));
end