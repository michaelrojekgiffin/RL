function opponent_parameters = empirical_priors_mg_simVsub_05_09_2017(OFFERS, predprey)
% this function estimates the true acceptance threshold of the distribution
% against which the opponent is playing. This is a logistic funciton
% consisting of two parameters, an intercept (paraneters(1)) and a slope
% (parameters(2)), both of which are estimated with this function and
% output in opponent parameters. 
% This function is necessary in order to simulate data using the parameters
% estimated from real subjects. The purpose of this fucntion is to estimate
% the true acceptance threshold of the distribution against which the
% subject was playing, so that we can write the simulation to be attempting
% to learn the same acceptance function, and then we can compare the
% behavior of our subject with the simulation playing with the same
% parameters trying to learn the same function, and hopefully they look the
% same.

ntr = 1; % nb of trials used per subject(1:first, 2:first two, 3:first three, etc...)

% compute frequences
sub_o = OFFERS(:);
sub_o_freq = zeros(1, 11);
for k = 0:10
    sub_o_freq(k+1)  = sum(sub_o==k);
end
sub_o_freq = sub_o_freq./(length(sub_o).*ntr);

% fit the parameters
options         = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);
n_rep           = 10;
parameters_rep  = NaN(n_rep,3);
ll_rep          = NaN(n_rep,1);
for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) PriorEstimation_mg_2017_09_04(x,sub_o, predprey),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
end
[~,pos]         = min(ll_rep);
parameters      = parameters_rep(pos(1),:);
ll              = ll_rep(pos(1),:);

opponent_parameters = parameters(1:2);


