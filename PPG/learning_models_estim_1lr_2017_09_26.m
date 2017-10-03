% This function generates the likelihood of each model/paramters.
% IMPORTANTLY, this is only meant to be run for models 1 and 3, since they
% have 1 learning rate. ALSO, this script appears to now be obsolete
% because it looks like I fixed the bug that prevented the other more
% general script from failing. I tried to make a generalizeable function
% that can have a variable number of free parameters but if the function
% ever has the option of having more parameters, even if I only assign
% those parameters values when certain criteria is met which isn't met in
% case of exeuction (e.g. model ~= 1), the minimal -loglik is still
% accomplished by estimating values for all free parameters and therefore
% I'm having separate functions for the different number of free parameters

function negLL = learning_models_estim_1lr_2017_09_26(params,o,r,a0,b0,nmodel, predprey)
% comp params

beta1 = params(1); % choice temperature
lr1   = params(2); % supraliminal learning rate
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
    a_t     = NaN(ntrial,1);             % logit intercept, updated at each trial
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
         % EV is different for predator and prey
        switch predprey
            case 'prey'
                EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);                       % compute EV of the offers given current model
            case 'predator'
                EV(ktrial,:)     = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
        end
        
        pc                  = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        lik(ktrial,kcond)   = pc(o(ktrial,kcond)+1);
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)      = r(ktrial) - EPc(ktrial);
                a_t(ktrial+1)   = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
%             case 2
%                 PE(ktrial)       = r(ktrial) - EPc(ktrial);
%                 if PE(ktrial)>0
%                     a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);
%                 else
%                     a_t(ktrial+1)  = a_t(ktrial) + 10.*lr2.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
%                 end
            case 3
                PE(ktrial)       = (r(ktrial) - EPc(ktrial)).*o(ktrial,kcond);
                a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
%             case 4
%                 PE(ktrial)       = (r(ktrial) - EPc(ktrial)).*o(ktrial,kcond);
%                 if PE(ktrial)>0
%                     a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);
%                 else
%                     a_t(ktrial+1)  = a_t(ktrial) + lr2.*PE(ktrial) ;                % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
%                 end
        end
    end
end

negLL = - sum(log(lik(:)));

% lik = -lik;                                                                               % loglikelyhood vector