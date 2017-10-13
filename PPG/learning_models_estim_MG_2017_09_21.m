% This function generates the likelihood of each model/paramters

function [negLL, EV, PA, a_t, pc, PE] = learning_models_estim_MG_2017_09_21(params,o,r,a0,b0,nmodel, predprey, R_o)
% comp params

switch nmodel
    case 1
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
%         lr2   = params(3); % supraliminal learning rate
        
    case 2
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
        
    case 3
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
%         lr2   = params(3); % supraliminal learning rate
        
    case 4
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
end
% beta1 = params(1); % choice temperature
% lr1   = params(2); % supraliminal learning rate
% lr2   = params(3); % supraliminal learning rate

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ntrial = size(o,1);
ncond  = size(o,2);

lik = NaN(ntrial,ncond);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));



for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t      = NaN(ntrial,1);            % logit intercept, updated at each trial
    pc      = NaN(ntrial, (numel(offers))); % expected offer, probability of offer being made
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
         % EV is different for predator and prey
        switch predprey
            case 'prey'
                EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);                       % compute EV of the offers given current model
                reward              = (10-o(ktrial,kcond));                           % for models 3 and 4
            case 'predator'
                EV(ktrial,:)     = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
                reward           = 10 - R_o(ktrial, kcond);
                % I think the following if statement  part should not be
                % used because I think that since the predator's left over
                % endowment is assured, it does not need to be included in
                % the calculation of prediction since it's not possible for
                % there to be any prediction error on this sum
%                 if r(ktrial) > 0
%                     reward          = (10-o(ktrial,kcond)) + (10 - R_o(ktrial, kcond));      % for use in models 3 and 4
%                 else
%                     reward          = (10-o(ktrial,kcond));
%                 end
        end
        
        pc(ktrial, :)       = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        lik(ktrial,kcond)   = pc(ktrial, o(ktrial,kcond)+1);
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)       = r(ktrial) - EPc(ktrial);
                a_t(ktrial+1)   = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 2
                PE(ktrial)       = r(ktrial) - EPc(ktrial);
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr2.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
            case 3
                PE(ktrial)       = (r(ktrial) - EPc(ktrial)).* reward;
                a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 4
                PE(ktrial)       = (r(ktrial) - EPc(ktrial)).*reward;
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + lr2.*PE(ktrial) ;                % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
        end
    end
end

negLL = - sum(log(lik(:)));

% lik = -lik;                                                                               % loglikelyhood vector
