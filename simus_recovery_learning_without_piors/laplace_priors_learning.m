% This function provide the computational model with paramters priors
function [post]=laplace_priors_learning(params,o,r,nmodel)

% log prior of parameters
beta1   = params(1);
a0      = params(2);
b0      = params(3);
lr1     = params(4);                                                            % choice temperature
lr2     = params(5);  

% [1 .5 .5 .5 0]
% pbeta1  = log(gampdf(beta1,1.2,5.0)); 
pbeta1  = log(gampdf(beta1,1.2,8.0)); 
pa0     = log(normpdf(a0,-6,2));
pb0     = log(betapdf(b0/5,1.1,1.1)); 
plr1    = log(betapdf(lr1,1.1,1.1));
plr2    = log(betapdf(lr2,1.1,1.1));                                        % the parameters are distrubution with mean + variance (different shapes) (?)


%%.5 .5 -.5
if nmodel==1
p = [pbeta1 pa0 pb0 plr1];
elseif nmodel==2
p = [pbeta1 pa0 pb0 plr1 plr2];
elseif nmodel==3
p = [pbeta1 pa0 pb0 plr1];
elseif nmodel==4
p = [pbeta1 pa0 pb0 plr1 plr2];
end

p = -sum(p);

l = learning_models_estim(params,o,r,nmodel);

post = p + l;


