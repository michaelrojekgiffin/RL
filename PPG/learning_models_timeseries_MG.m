% This function generates the likelihood of each model/paramters

function [o,r,pe,at,pd, R_o] = learning_models_timeseries_MG(paramsP,paramsR,ntrial,a0,b0,nmodel, predprey)
% comp params
beta1 = paramsP(1); % choice temperature
lr1   = paramsP(2); % supraliminal learning rate
lr2   = paramsP(3); % supraliminal learning rate

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

R_o = NaN(ntrial,ncond);               % the actual offers of the rival, given the model
%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial+1,1);           % logit intercept, updated at each trial
    
    % pre-allocate for rival, this is to calculate reward prediction errors
    %in predators
    R_PA    = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    R_EV    = NaN(ntrial,numel(offers)); % expected value
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
         PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        
        %         R_PA(ktrial,:)  = logitp([Ra(kcond),Rb(kcond)],offers);    % Do the same for the rival (opponent)
        % the above is commented out because in this case all the different
        % conditions will be learning the same distribution
        R_PA(ktrial,:)  = logitp([Ra,Rb],offers);    % Do the same for the rival (opponent)
        
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
                %                 R_EV(ktrial,:)   = (endow - offers) +((endow - offers).*R_PA(ktrial,:));  % compute rivals EV
                
                tempEV = R_EV(ktrial,:); % predator EV is calculated a bit different
                for tc = 1:length(tempEV)
                    tempEV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc))); 
                end
                R_EV(ktrial,:)   = tempEV;
                
            case 'predator'
%                 EV(ktrial,:)     = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
                R_EV(ktrial,:)   = (endow - offers).*R_PA(ktrial,:);                     % compute rivals EV
                
                tempEV = EV(ktrial,:);
                for tc = 1:length(tempEV)
                    tempEV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc))); 
                end
                tempEV(1) = 10;
                
%                 EV(ktrial,:)        = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
                EV(ktrial, :)       = tempEV;
               
        end
        
        pc                  = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        pd                  = makedist('multinomial','probabilities',pc);                  % estimate the pdf from the pC
        kO                  = random(pd);                                                  % selected offer, in the 1:numel(offer) spavce
        
        o(ktrial,kcond)     = offers(kO);                                                  % resample Offer in pdf (="soft-max")
        
        Pd                  = logitp([Ra(kcond),Rb(kcond)],o(ktrial,kcond));               % Reciever Estimated accepantce proba of accepting the offer
        
        
        R_pc                = exp(Rb(kcond).*R_EV(ktrial,:)) ./ sum(exp(Rb(kcond).*R_EV(ktrial,:)));   % multinomial choice function of rival
        R_pd                = makedist('multinomial','probabilities',R_pc);                            % estimate the pdf from the R_pC
        R_kO                = random(R_pd);                                                            % selected offer, in the 1:numel(offer) spavce
        R_o(ktrial,kcond)   = offers(R_kO);    
        
        
        % for predator offer must be less than opponent's offer, for prey
        % offer must be less than or equal to opponent's offer
        switch predprey 
            case 'predator'
                r(ktrial,kcond)     = double(rand(1)<Pd);                             % Sampling Reciever's decision given the proba.
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
%                 
            case 'prey'
                r(ktrial,kcond)     = double(rand(1)<=Pd);                            % Sampling Reciever's decision given the proba.
                reward              = (10-o(ktrial,kcond));                           % for models 3 and 4
        end
%         r(ktrial,kcond)     = double(rand(1)<Pd);                                   % Sampling Reciever's decision given the proba.
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)         = r(ktrial,kcond) - EPc(ktrial);
                a_t(ktrial+1)      = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   
                
            case 2
                PE(ktrial)         = r(ktrial,kcond) - EPc(ktrial);
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr2.*PE(ktrial);                  
                end
            case 3
                PE(ktrial)         = (r(ktrial,kcond) - EPc(ktrial)).*reward;
                a_t(ktrial+1)      = a_t(ktrial) + lr1.*PE(ktrial);                   
                
            case 4
                PE(ktrial)         = (r(ktrial,kcond) - EPc(ktrial)).*reward;
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + lr2.*PE(ktrial);                          
                end
        end
        
    end
    
    
    pe(:,kcond) = PE;
    at(:,kcond) = a_t; % estimate of opponent's choice function slope on each trial
end


% lik = -lik;                                                                               % loglikelyhood vector
