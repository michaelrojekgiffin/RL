% This function provide the computational model with paramters priors
function [post]=laplace_priors_priors(params,tr, o)
% tr is the array of the binary decisions, and o if the offer selected, I
% think it's 0 for offer 1, and 1 for offer 2, not totally sure though

% log prior of parameters

a1    = params(1);                                        % choice temperature
b1    = params(2);                                        % factual learning rate
beta1 = params(3);

% these are constrained within values based on the literature, ask Mael for
% more details about how he arrived at these
% [1 .5 .5 .5 0]
pa1    = log(normpdf(a1, -5,2));
pb1    = log(gampdf(b1,1.2,5.0));                         % the parameters are distrubution with mean + variance (different shapes) (?)
pbeta1 = log(gampdf(beta1,1.2,5.0));                      % the parameters are distrubution with mean + variance (different shapes) (?)


% plr1   = log(betapdf(lr1, 1.1,1.1));
%%.5 .5 -.5

p = [pa1 pb1 pbeta1];
p = -sum(p);

l = PriorEstimation(params,tr, o);

post = p + l;


