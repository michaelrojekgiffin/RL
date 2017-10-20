% This function generates the likelihood of each model/paramters
% this model was adapted in order to use in the model_vs_subject script, so
% that I can use the actual opponent distribution instead of simulating it.

function [o_all,r_all,pe_all,at_all,pc_all, R_o_all, PA_all, EV_all] = learning_models_timeseries_MG_2017_10_03(paramsP,paramsR,ntrial,a0,b0,nmodel, predprey, opponent_o)
% comp params
beta1 = paramsP(1); % choice temperature
lr1   = paramsP(2); % supraliminal learning rate
lr2   = paramsP(3); % supraliminal learning rate

Ra = paramsR(1);
Rb = paramsR(2);

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));

% this is a vestige from Mael's more generalizeable script, I'm only going
% to use this for one condition, so here we are
ncond  = size(opponent_o, 2); %size(paramsR,2);

o  = NaN(ntrial,ncond);
r  = NaN(ntrial,ncond);
pe  = NaN(ntrial,ncond);
at  = NaN(ntrial+1,ncond);

R_o = NaN(ntrial,ncond);               % the actual offers of the rival, given the model
%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% to save everything in one array
o_all       = NaN(ntrial*ncond,1);
r_all       = NaN(ntrial*ncond,1);
pe_all      = NaN(ntrial*ncond,1);
at_all      = NaN((ntrial*ncond)+1,1);
pc_all      = NaN(ntrial*ncond, (numel(offers))); % expected offer, probability of offer being made
R_o_all     = NaN(ntrial*ncond,1);               % the actual offers of the rival, given the model
PA_all      = NaN(ntrial*ncond,numel(offers)); % estimated probability of accepting all offers
EV_all      = NaN(ntrial*ncond,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
pc_counter  = 1;

for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial+1,1);           % logit intercept, updated at each trial
    EVo     = NaN(ntrial, 1);            % EV of the selected offer
    
    % pre-allocate for rival, this is to calculate reward prediction errors
    %in predators
    R_PA    = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    R_EV    = NaN(ntrial,numel(offers)); % expected value
    
    pc      = NaN(ntrial, (numel(offers))); % expected offer, probability of offer being made
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
        
        %         R_PA(ktrial,:)  = logitp([Ra(kcond),Rb(kcond)],offers);    % Do the same for the rival (opponent)
        % the above is commented out because in this case all the different
        % conditions will be learning the same distribution
        R_PA(ktrial,:)  = logitp([Ra,Rb],offers);    % Do the same for the rival (opponent)
        % EV is different for predator and prey
        switch predprey
            case 'prey'
                EV(ktrial,:)     = (endow - offers).* PA(ktrial,:);                       % compute EV of the offers given current model
                R_EV(ktrial,:)   = (endow - offers) +((endow - offers).*R_PA(ktrial,:));  % compute rivals EV
            case 'predator'
                EV(ktrial,:)     = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
                R_EV(ktrial,:)   = (endow - offers).*R_PA(ktrial,:);                     % compute rivals EV
        end
        
        pc(ktrial, :)        = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        pd                  = makedist('multinomial','probabilities',pc(ktrial, :));                  % estimate the pdf from the pC
        kO                  = random(pd);                                                  % selected offer, in the 1:numel(offer) spavce
        
        o(ktrial,kcond)     = offers(kO);                                                  % resample Offer in pdf (="soft-max")
        
%         Pd                  = logitp([Ra(kcond),Rb(kcond)],o(ktrial,kcond));               % Reciever Estimated accepantce proba of accepting the offer
        Pd                  = logitp([Ra,Rb],o(ktrial,kcond));               % Reciever Estimated accepantce proba of accepting the offer
       
        
%         R_pc                = exp(Rb(kcond).*R_EV(ktrial,:)) ./ sum(exp(Rb(kcond).*R_EV(ktrial,:)));   % multinomial choice function of rival
        R_pc                = exp(Rb.*R_EV(ktrial,:)) ./ sum(exp(Rb.*R_EV(ktrial,:)));   % multinomial choice function of rival
        
        R_pd                = makedist('multinomial','probabilities',R_pc);                            % estimate the pdf from the R_pC
        R_kO                = random(R_pd);                                                            % selected offer, in the 1:numel(offer) spavce
        R_o(ktrial,kcond)   = offers(R_kO);    
        
        
        % for predator offer must be less than opponent's offer, for prey
        % offer must be less than or equal to opponent's offer
        switch predprey 
            case 'predator'
                r(ktrial,kcond)     = double(rand(1)<Pd);                             % Sampling Reciever's decision given the proba.
                reward           = 10 - opponent_o(ktrial, kcond);
%                 reward           = 10 - R_o(ktrial, kcond);
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
    
    % storing everything in a single array
    o_all(pc_counter:pc_counter+ntrial-1, :)    = o(:, kcond);
    r_all(pc_counter:pc_counter+ntrial-1, :)    = r(:, kcond);
    pe_all(pc_counter:pc_counter+ntrial-1, :)   = pe(:, kcond);
    
    pc_all(pc_counter:pc_counter+ntrial-1, :)   = pc; 
    R_o_all(pc_counter:pc_counter+ntrial-1, :)  = R_o(:, kcond);
    PA_all(pc_counter:pc_counter+ntrial-1, :)   = PA;
    EV_all(pc_counter:pc_counter+ntrial-1, :)   = EV;
    
    if pc_counter == 1
        at_all(pc_counter:pc_counter+ntrial, :)         = at(:, kcond); % has an extra row due to the prior
    else
        at_all(pc_counter+1:pc_counter+ntrial, :)     = at(2:end, kcond); % has an extra row due to the prior
    end
    
    pc_counter = pc_counter + ntrial;
end


% lik = -lik;                                                                               % loglikelyhood vector
