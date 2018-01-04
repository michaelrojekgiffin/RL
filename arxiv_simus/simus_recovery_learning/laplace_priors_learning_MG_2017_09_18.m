% This function provide the computational model with paramters priors
function [post]=laplace_priors_learning_MG_2017_09_18(params,o,r,a0,b0,nmodel)

% log prior of parameters
beta1  = params(1);
lr1    = params(2);                                                            % choice temperature
if length(params) > 2
    lr2    = params(3);
end

% [1 .5 .5 .5 0]
% pbeta1  = log(gampdf(beta1,1.2,5.0)); 
pbeta1  = log(gampdf(beta1,1.2,8.0)); 
plr1    = log(betapdf(lr1,1.1,1.1));
if length(params) > 2
    plr2    = log(betapdf(lr2,1.1,1.1));                                        % the parameters are distrubution with mean + variance (different shapes) (?)
end

%%.5 .5 -.5
if nmodel==1
    p = [pbeta1 plr1];
    l = learning_models_estim_1lr_MG_2017_09_18(params,o,r,a0,b0,nmodel);
elseif nmodel==2
    p = [pbeta1 plr1 plr2];
    l = learning_models_estim_2lr_MG_2017_09_18(params,o,r,a0,b0,nmodel);
elseif nmodel==3
    p = [pbeta1 plr1];
    l = learning_models_estim_1lr_MG_2017_09_18(params,o,r,a0,b0,nmodel);
elseif nmodel==4
    p = [pbeta1 plr1 plr2];
    l = learning_models_estim_2lr_MG_2017_09_18(params,o,r,a0,b0,nmodel);
end

p = -sum(p);

% l = learning_models_estim(params,o,r,a0,b0,nmodel);

post = p + l;


