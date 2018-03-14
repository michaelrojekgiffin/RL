% This function generates the likelihood of each model/paramters

function [o,r, opponent_o, pe,at,pd] = learning_models_timeseries(paramsP,paramsR,ntrial,nmodel, predprey)
% comp params
beta1 = paramsP(1); % choice temperature
a0    = paramsP(2);
b0    = paramsP(3);
lr1   = paramsP(4); % supraliminal learning rate
lr2   = paramsP(5); % supraliminal learning rate

Ra = paramsR(1,:);
Rb = paramsR(2,:);

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ncond  = size(paramsR,2);

o  = NaN(ntrial,ncond);
r  = NaN(ntrial,ncond);
pe  = NaN(ntrial,ncond);
at  = NaN(ntrial+1,ncond);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial+1,1);            % logit intercept, updated at each trial
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        R_PA(ktrial,:)  = logitp([Ra(kcond),Rb(kcond)],offers);                       % Do the same for the rival (opponent)
        
        % importantly, PA represents the probability that the opponent
        % chooses a lower investment than the subject. In order to
        % calculate EV for predator, I need a probability value for each investment.
        
        rPA     = PA(ktrial, :);
        for rrr = 1:length(rPA)
            rPA(rrr)     = rPA(rrr) - sum(rPA(1:rrr-1));
        end
        
        % importantly, PA represents the probability that the opponent
        % chooses a lower investment than the subject. In order to
        % calculate EV for predator, I need a probability value for each investment.
        rPA     = PA(ktrial, :);
        for rrr = 1:length(rPA)
            rPA(rrr)     = rPA(rrr) - sum(rPA(1:rrr-1));
        end
        % EV is different for predator and prey
        switch predprey
            case 'prey'
                EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);                       % compute EV of the offers given current model

%                 tempEV = R_EV(ktrial,:); % predator EV is calculated a bit different
%                 for tc = 1:length(tempEV)
%                     tempEV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc))); 
%                 end
%                 R_EV(ktrial,:)   = tempEV;
                
            case 'predator'
                
                R_EV(ktrial,:)   = (endow - offers).*R_PA(ktrial,:);                     % compute rivals EV
                
                tempEV = EV(ktrial,:);
                for tc = 1:length(tempEV)
                    tempEV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc))); 
                end
                tempEV(1) = 10;
                
                EV(ktrial, :)       = tempEV;
               
        end
        
        pc      = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        pd      = makedist('multinomial','probabilities',pc);           % estimate the pdf from the pC
        kO      = random(pd);                                           % selected offer, in the 1:numel(offer) spavce
        
        o(ktrial,kcond)     = offers(kO);                                % resample Offer in pdf (="soft-max")
        
        Pd                  = logitp([Ra(kcond),Rb(kcond)],o(ktrial,kcond));    % Reciever Estimated accepantce proba of accepting the offer
        
        % for predator offer must be less than opponent's offer, for prey
        % offer must be less than or equal to opponent's offer
        switch predprey 
            case 'predator'
                % opponent_o is the opponent's offer, which I'm simulating
                % here based on the actual probabilities of the
                % computerized opponent prey in the hor dataset.
                opponent_o(ktrial,kcond)     = randsample([0:10],1,true, [0.066666667 0.064509804 0.074901961 0.089215686 0.109411765 0.181568627 0.161764706 0.146274510 0.077647059 0.021372549 0.006666667]);
                r(ktrial,kcond)              = double(rand(1)<Pd);                             % Sampling Reciever's decision given the proba.
                reward                       = 10 - opponent_o(ktrial,kcond);

            case 'prey'
                % opponent_o is the opponent's offer, which I'm simulating
                % here based on the actual probabilities of the
                % computerized opponent predator in the hor dataset.
                opponent_o(ktrial,kcond)     = randsample([0:10],1,true, [0.302880658 0.095884774 0.058024691 0.085185185 0.105967078 0.116255144 0.148148148 0.033950617 0.033744856 0.014814815 0.005144033]);
                r(ktrial,kcond)     = double(rand(1)<=Pd);                            % Sampling Reciever's decision given the proba.
                reward              = (10-o(ktrial,kcond));                           % for models 3 and 4
        end
        
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)       = r(ktrial,kcond) - EPc(ktrial);
                a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 2
                PE(ktrial)       = r(ktrial,kcond) - EPc(ktrial);
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr2.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
            case 3
                PE(ktrial)       = (r(ktrial,kcond) - EPc(ktrial)).*reward;
                a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 4
                PE(ktrial)       = (r(ktrial,kcond) - EPc(ktrial)).*reward;
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + lr2.*PE(ktrial);                          % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
        end
        
    end
    
    
    pe(:,kcond) = PE;
    at(:,kcond) = a_t;
end


% lik = -lik;                                                                               % loglikelyhood vector
