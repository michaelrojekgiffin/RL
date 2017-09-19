% This function generates the likelihood of each model/paramters

function negLL = learning_model_2params(params,o,r,a0,b0)
% comp params
beta1 = params(1); % choice temperature
lr1   = params(2); % supraliminal learning rate

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ntrial = size(o,1);
ncond  = size(o,2);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));



negLL = 0;

for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    CPE     = NaN(ntrial,1);             % Choice prediction error
    a_t      = NaN(ntrial,1);            % logit intercept, updated at each trial
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);               % compute EV of the offers given current model
        
        pc = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        
        lik = pc(o(ktrial,kcond)+1);
        
        EPc(ktrial)      = logitp([a_t(ktrial),b0],o(ktrial,kcond));
        CPE(ktrial)      = r(ktrial) - EPc(ktrial);
        a_t(ktrial+1,:)  = a_t(ktrial) + lr1.*CPE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        
        negLL = negLL - log(lik);
    end
end

% lik = -lik;                                                                               % loglikelyhood vector
