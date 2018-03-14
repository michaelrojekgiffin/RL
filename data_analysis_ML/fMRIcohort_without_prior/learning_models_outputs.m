% This function generates the likelihood of each model/paramters

function [PE,eO,a_t] = learning_models_outputs(params,o,r,nmodel)
% comp params

beta1   = params(1); % choice temperature
a0      = params(2); % supraliminal learning rate
b0      = params(3); % supraliminal learning rate
lr1     = params(4); % supraliminal learning rate
lr2     = params(5); % supraliminal learning rate


% task param
offers  = 0:1:20;
endow   = 20*ones(1,numel(offers));
ntrial = size(o,1);
ncond  = size(o,2);

lik = NaN(ntrial,ncond);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


PE      = NaN(ntrial,ncond);            % Choice prediction error
a_t     = NaN(ntrial,ncond);            % logit intercept, updated at each trial
eO      = NaN(ntrial,ncond);            % expected value of the offer

for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    
    pc      = NaN(ntrial, (numel(offers))); % expected offer, probability of offer being made
    
    %initiate
    a_t(1,kcond)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)        = logitp([a_t(ktrial,kcond),b0],offers);                           % compute proba of accepting the offers given current model
        EV(ktrial,:)        = (endow - offers).* PA(ktrial,:);                                 % compute EV of the offers given current model        
        pc(ktrial,:)        = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));    % multinomial choice function
        eO(ktrial,kcond)    = pc(ktrial,:)*offers';
        lik(ktrial,kcond)   = pc(ktrial,o(ktrial,kcond)+1);
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                % PE(ktrial)       = r(ktrial) - EPc(ktrial);
                PE(ktrial,kcond)       = r(ktrial,kcond) - EPc(ktrial);
                a_t(ktrial+1,kcond)   = a_t(ktrial,kcond) + 20.*lr1.*PE(ktrial,kcond);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 2
                % PE(ktrial)       = r(ktrial) - EPc(ktrial);
                PE(ktrial,kcond)       = r(ktrial,kcond) - EPc(ktrial);
                if PE(ktrial,kcond)>0
                    a_t(ktrial+1,kcond)  = a_t(ktrial,kcond) + 20.*lr1.*PE(ktrial,kcond);
                else
                    a_t(ktrial+1,kcond)  = a_t(ktrial,kcond) + 20.*lr2.*PE(ktrial,kcond);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
            case 3
                % PE(ktrial)       = (r(ktrial) - EPc(ktrial)).*o(ktrial,kcond);
                PE(ktrial,kcond)       = (r(ktrial,kcond) - EPc(ktrial)).*(20-o(ktrial,kcond));
                a_t(ktrial+1,kcond)  = a_t(ktrial,kcond) + lr1.*PE(ktrial,kcond);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 4
                % PE(ktrial)       = (r(ktrial) - EPc(ktrial)).*o(ktrial,kcond);
                PE(ktrial,kcond)       = (r(ktrial,kcond) - EPc(ktrial)).*(20-o(ktrial,kcond));
                if PE(ktrial,kcond)>0
                    a_t(ktrial+1,kcond)  = a_t(ktrial,kcond) + lr1.*PE(ktrial,kcond);
                else
                    a_t(ktrial+1,kcond)  = a_t(ktrial,kcond) + lr2.*PE(ktrial,kcond) ;                % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
        end
    end
end
