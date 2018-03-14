% This function generates the likelihood of each model/paramters

function negLL = learning_models_estim(params,o,r,nmodel)
% comp params

beta1nS = params(1); % choice temperature
a0nS    = params(2); % supraliminal learning rate
b0nS    = params(3); % supraliminal learning rate
lr1nS   = params(4); % supraliminal learning rate
lr2nS   = params(5); % supraliminal learning rate

switch nmodel
    case 1
    case 2
        beta1S = params(6);
    case 3
        a0S = params(6);
    case 4
        b0S = params(6);
    case 5
        lr1S = params(6);
    case 6
        lr2S = params(6);
    case 7
        beta1S = params(6);
        a0S = params(7);
    case 8
        beta1S = params(6);
        b0S = params(7);
    case 9
        beta1S = params(6);
        lr1S = params(7);
    case 10
        beta1S = params(6);
        lr2S = params(7);
    case 11
        a0S = params(6);
        b0S = params(7);
    case 12
        a0S = params(6);
        lr1S = params(7);
    case 13
        a0S = params(6);
        lr2S = params(7);
    case 14
        b0S = params(6);
        lr1S = params(7);
    case 15
        b0S = params(6);
        lr2S = params(7);
    case 16
        lr1S = params(6);
        lr2S = params(7);
    case 17
        beta1S = params(6);
        a0S = params(7);
        b0S = params(8);
    case 18
        beta1S = params(6);
        a0S = params(7);
        lr1S = params(8);
    case 19
        beta1S = params(6);
        a0S = params(7);
        lr2S = params(8);
    case 20
        beta1S = params(6);
        b0S = params(7);
        lr1S = params(8);
    case 21
        beta1S = params(6);
        b0S = params(7);
        lr2S = params(8);
    case 22
        beta1S = params(6);
        lr1S = params(7);
        lr2S = params(8);
    case 23
        a0S = params(6);
        b0S = params(7);
        lr1S = params(8);
    case 24
        a0S = params(6);
        b0S = params(7);
        lr2S = params(8);
    case 25
        a0S = params(6);
        lr1S = params(7);
        lr2S = params(8);
    case 26
        b0S = params(6);
        lr1S = params(7);
        lr2S = params(8);
    case 27
        beta1S = params(6); % choice temperature
        a0S    = params(7); % supraliminal learning rate
        b0S    = params(8); % supraliminal learning rate
        lr1S   = params(9); % supraliminal learning rate
    case 28
        beta1S = params(6); % choice temperature
        a0S    = params(7); % supraliminal learning rate
        b0S    = params(8); % supraliminal learning rate
        lr2S   = params(9); % supraliminal learning rate
    case 29
        beta1S = params(6); % choice temperature
        a0S    = params(7); % supraliminal learning rate
        lr1S   = params(8); % supraliminal learning rate
        lr2S   = params(9); % supraliminal learning rate
    case 30
        beta1S = params(6); % choice temperature
        b0S    = params(7); % supraliminal learning rate
        lr1S   = params(8); % supraliminal learning rate
        lr2S   = params(9); % supraliminal learning rate
    case 31
        a0S    = params(6); % supraliminal learning rate
        b0S    = params(7); % supraliminal learning rate
        lr1S   = params(8); % supraliminal learning rate
        lr2S   = params(9); % supraliminal learning rate
    case 32
        beta1S = params(6); % choice temperature
        a0S    = params(7); % supraliminal learning rate
        b0S    = params(8); % supraliminal learning rate
        lr1S   = params(9); % supraliminal learning rate
        lr2S   = params(10); % supraliminal learning rate
end


% task param
offers  = 0:1:20;
endow   = 20*ones(1,numel(offers));
ntrial = size(o,1);
ncond  = size(o,2);

lik = NaN(ntrial,ncond);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));





for kcond = 1:ncond
    S = double(kcond>6);
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial,1);            % logit intercept, updated at each trial
    pc      = NaN(ntrial, (numel(offers))); % expected offer, probability of offer being made
    
    
    beta1 = beta1nS;
    a0 = a0nS;
    b0 = b0nS;
    lr1 = lr1nS;
    lr2 = lr2nS;
    
    switch nmodel
        
        case 1
        case 2
            beta1 = (1-S)*beta1nS + S*(beta1S);
        case 3
            a0 = (1-S)*a0nS + S*(a0S);
        case 4
            b0 = (1-S)*b0nS + S*(b0S);
        case 5
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 6
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 7
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
        case 8
            beta1 = (1-S)*beta1nS + S*(beta1S);
            b0 = (1-S)*b0nS + S*(b0S);
        case 9
            beta1 = (1-S)*beta1nS + S*(beta1S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 10
            beta1 = (1-S)*beta1nS + S*(beta1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 11
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
        case 12
            a0 = (1-S)*a0nS + S*(a0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 13
            a0 = (1-S)*a0nS + S*(a0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 14
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 15
            b0 = (1-S)*b0nS + S*(b0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 16
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 17
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
        case 18
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 19
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 20
            beta1 = (1-S)*beta1nS + S*(beta1S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 21
            beta1 = (1-S)*beta1nS + S*(beta1S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 22
            beta1 = (1-S)*beta1nS + S*(beta1S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 23
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 24
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 25
            a0 = (1-S)*a0nS + S*(a0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 26
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 27
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
        case 28
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 29
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 30
            beta1 = (1-S)*beta1nS + S*(beta1S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 31
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
        case 32
            beta1 = (1-S)*beta1nS + S*(beta1S);
            a0 = (1-S)*a0nS + S*(a0S);
            b0 = (1-S)*b0nS + S*(b0S);
            lr1 = (1-S)*lr1nS + S*(lr1S);
            lr2 = (1-S)*lr2nS + S*(lr2S);
    end
    
    
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        PA(ktrial,:)        = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        EV(ktrial,:)        = (endow - offers).* PA(ktrial,:);               % compute EV of the offers given current model
        pc(ktrial,:)        = exp(beta1*EV(ktrial,:)) ./ sum(exp(beta1*EV(ktrial,:)));   % multinomial choice function
        lik(ktrial,kcond)   = pc(ktrial,o(ktrial,kcond)+1);
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        PE(ktrial)          = r(ktrial,kcond) - EPc(ktrial);
        PEval               = PE(ktrial)>0;
        a_t(ktrial+1)       = a_t(ktrial) + 20*(lr1*PEval + lr2*(~PEval))*PE(ktrial);
    end
    
end

negLL = - sum(log(lik(:)));

% lik = -lik;                                                                               % loglikelyhood vector
