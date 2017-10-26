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

rPA     = PA;
for rrr = 1:length(rPA)
    rPA(rrr)     = rPA(rrr) - sum(rPA(1:rrr-1));
end

if strcmp(predprey, 'prey')
    EV      = (endow - offers).* PA;           % compute EV of the offers given current model
else
    %     EV      = (endow - offers) + (endow - offers).* PA;           % compute EV of the offers given current model
    EV = NaN(1, length(offers));
    for tc = 1:length(offers)
        EV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc)));
    end
    
    %                 EV(ktrial,:)        = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
%     EV(ktrial, :)       = tempEV;
    
    
end
pc      = exp(betaX.*EV) ./ sum(exp(betaX.*EV));

%likelihood
negLL =  sum(-log(pc(o+1)));