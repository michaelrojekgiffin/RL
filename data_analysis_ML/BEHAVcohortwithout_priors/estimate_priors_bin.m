% This function generates the likelihood of each model/paramters

function negLL = estimate_priors_bin(params,tr,o)
% comp params
beta0 = params(1);  % intercept of estimate                     % choice temperature
beta1 = params(2);  % slope of estimate
betaX = params(3);  % temperature of logit (function of beta0 & beta1) update 

% task param
offers  = 0:1:20;
endow   = 20*ones(1,numel(offers));

pc = NaN(numel(o),1);

%funct
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Compute the p(offers)
PA      = logitp([beta0,beta1],offers);        % compute proba of accepting the offers given current model
EV      = (endow - offers).* PA;               % compute EV of the offers given current model

for k = 1:length(o)
    if ~isnan(o(k))
    pc(k)     = exp(betaX.*EV(tr(k,1)+1) ) ./ sum(exp(betaX.*EV(tr(k,:)+1)));   
    end % resample Offer in pdf (="soft-max")   
end

nfin = sum(~isnan(o));
ll(1:nfin) = o(1:nfin) .*pc(1:nfin) + (1-o(1:nfin) ).*(1-pc(1:nfin));

tol = 1e-12;
ll(ll<tol) = tol;
ll(ll>1-tol) = 1-tol;

%likelihood
negLL =  sum(-log(ll));