function eval = evaluate_step(priorL,newL)
% evaluates whether a step should be taken
if newL >= priorL
    eval = 1;
else
    m = rand(1);
    p = exp(newL-priorL);
    if m >= p
        eval = 0;
    else
        eval = 1;
    end
end

end