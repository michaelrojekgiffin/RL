% This function generates the likelihood of each model/paramters

function negLL = PriorEstimation_mg_2017_09_04(params,o, predprey)
% comp params
beta0 = params(1);  % intercept of estimate                     % choice temperature
beta1 = params(2);  % slope of estimate
betaX = params(3);  % temperature of logit (function of beta0 & beta1) update 

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));

%funct
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Compute the p(offers)
PA      = logitp([beta0,beta1],offers);        % compute proba of accepting the offers given current model
if strcmp(predprey, 'prey')
    EV      = (endow - offers).* PA;           % compute EV of the offers given current model
else
    EV      = (endow - offers) + (endow - offers).* PA;           % compute EV of the offers given current model
end
pc      = exp(betaX.*EV) ./ sum(exp(betaX.*EV));

%likelihood
negLL =  sum(-log(pc(o+1)));