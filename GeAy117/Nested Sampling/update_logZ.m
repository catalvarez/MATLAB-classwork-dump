function new_logZ = update_logZ(current_logZ,log_Ak)
% adds new area to Z in log-space

if current_logZ >= log_Ak
    new_logZ = current_logZ + log(1+exp(log_Ak-current_logZ));
else
    new_logZ = log_Ak+log(1+exp(current_logZ-log_Ak));
end
end