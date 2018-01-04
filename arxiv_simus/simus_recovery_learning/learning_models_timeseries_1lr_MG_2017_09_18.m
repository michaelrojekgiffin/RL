% This function generates the likelihood of each model/paramters

function [o,r] = learning_models_timeseries_1lr_MG_2017_09_18(paramsP,paramsR,ntrial,a0,b0,nmodel)
% comp params
beta1 = paramsP(1); % choice temperature
lr1   = paramsP(2); % supraliminal learning rate
% lr2   = paramsP(3); % supraliminal learning rate

Ra = paramsR(1,:);
Rb = paramsR(2,:);

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ncond  = size(paramsR,2);

o = NaN(ntrial,ncond);
r = NaN(ntrial,ncond);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial,1);            % logit intercept, updated at each trial
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);             % compute EV of the offers given current model
        
        pc      = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        pd      = makedist('multinomial','probabilities',pc);           % estimate the pdf from the pC
        kO      = random(pd);                                           % selected offer, in the 1:numel(offer) spavce
        
        o(ktrial,kcond)     = offers(kO);                                % resample Offer in pdf (="soft-max")
        
        Pd                  = logitp([Ra(kcond),Rb(kcond)],o(ktrial,kcond));    % Reciever Estimated accepantce proba of accepting the offer
        r(ktrial,kcond)     = double(rand(1)<Pd);                               % Sampling Reciever's decision given the proba.
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)       = r(ktrial,kcond) - EPc(ktrial);
                a_t(ktrial+1,:)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 2
                PE(ktrial)       = r(ktrial,kcond) - EPc(ktrial);
                a_t(ktrial+1,:)  = a_t(ktrial) + 10.*lr1.*PE(ktrial).*(PE(ktrial)>0)+ 10.*lr2.*PE(ktrial).*(PE(ktrial)<0);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 3
                PE(ktrial)       = (r(ktrial,kcond) - EPc(ktrial)).*o(ktrial,kcond);
                a_t(ktrial+1,:)  = a_t(ktrial) + lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 4
                PE(ktrial)       = (r(ktrial,kcond) - EPc(ktrial)).*o(ktrial,kcond);
                a_t(ktrial+1,:)  = a_t(ktrial) + lr1.*PE(ktrial).*(PE(ktrial)>0)+ lr2.*PE(ktrial).*(PE(ktrial)<0);                          % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        end
        
    end
end


% lik = -lik;                                                                               % loglikelyhood vector
