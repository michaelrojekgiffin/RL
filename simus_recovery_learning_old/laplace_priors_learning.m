% This function provide the computational model with paramters priors
function [post]=laplace_priors_learning(params,o,r,a0,b0)

% log prior of parameters

a1    = params(1);                                                            % choice temperature
beta1 = params(2);

% [1 .5 .5 .5 0]
pa1    = log(normpdf(a1, -5,2));
pbeta1 = log(gampdf(beta1,1.2,5.0));                                            % the parameters are distrubution with mean + variance (different shapes) (?)


% plr1   = log(betapdf(lr1, 1.1,1.1));
%%.5 .5 -.5

p = [pa1 pbeta1];
p = -sum(p);

l = learning_model_2params(params,o,r,a0,b0);

post = p + l;


