function index = perturb_param(param)
% randomly determines which parameter to perturb, returns param and
% step-size
index = randi(length(param),1);
% i = ceil(rand(1)*length(param(1,:)));
% p = param(:,i);
% step = stepsize(i);
end