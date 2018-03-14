% This function generates the likelihood of each model/paramters

function [negLL, all_EV, all_PA, a_t, all_pc, all_PE, risk, risk_pe] = learning_models_estim_MG(params,o,r,a0,b0,nmodel, predprey, R_o)
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

% to keep everything in one array for the output
all_pc      = NaN(ncond*ntrial, (numel(offers))); 
all_PE      = NaN(ncond*ntrial, 1);
all_EV      = NaN(ncond*ntrial, (numel(offers)));
all_PA      = NaN(ncond*ntrial, (numel(offers)));
pc_counter  = 1;

risk        = NaN(ncond*ntrial, 1); 
risk_pe     = NaN(ncond*ntrial, 1); 
risk_count  = 0;

for kcond = 1:ncond
    % pre-allocate
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers% V       = NaN(ntrial); % estimated expected value of all offers
    EPc     = NaN(ntrial,1);             % estimated probability of accepting the offer
    PE      = NaN(ntrial,1);             % Choice prediction error
    a_t     = NaN(ntrial,1);            % logit intercept, updated at each trial
    pc      = NaN(ntrial, (numel(offers))); % expected offer, probability of offer being made
    
    %initiate
    a_t(1)     = a0;                    % initial value of teh logit intercept
    for ktrial = 1:ntrial
        risk_count = risk_count + 1;
        PA(ktrial,:)     = logitp([a_t(ktrial),b0],offers);            % compute proba of accepting the offers given current model
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
                EV(ktrial,:)        = (endow - offers).* PA(ktrial,:);                       % compute EV of the offers given current model
                reward              = (10-o(ktrial,kcond));                                  % for models 3 and 4
                % this is based on the formula Frans gave me, it sums over
                % the squares of all the possible prediction errors for the
                % subject's decision
                sspe                = (((10 - offers) - EV(ktrial,o(ktrial, kcond)+1)).^2);
                risk(risk_count)    = sum(rPA .* sspe);
                
%                 risk(risk_count)    = sum (rPA .* (((10 - offers) - EV(ktrial,o(ktrial, kcond)+1)).^2));
                risk_pe(risk_count) = ((10 - o(ktrial,kcond) - EV(ktrial, o(ktrial, kcond)+1)) ^2) - risk(risk_count) ;
            case 'predator'
                tempEV = EV(ktrial,:);
                for tc = 1:length(tempEV)
                    if tc == 1
                        tempEV(tc) = ((endow(tc) - offers(tc)));
                    else
                        tempEV(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (rPA(1:tc)));
                    end
                    
                end
%                 tempEV(1) = 10; % EV for an investment of 0 should always be 10
                
%                 EV(ktrial,:)        = (endow - offers) +((endow - offers).* PA(ktrial,:));   % compute EV of the offers given current model
                EV(ktrial, :)       = tempEV;
                
                reward              = 10 - R_o(ktrial, kcond);
                
                %                 if o(ktrial, kcond) == 0
                %                     risk(risk_count)    = 0;
                %                     risk_pe(risk_count) = 0;
                %                 else
                %                     sspe                = (((10 - offers) + 10 - o(ktrial,kcond)) - EV(ktrial, o(ktrial, kcond)+1)).^2;
                %                     risk(risk_count)    = sum(rPA .* sspe);
                %                     risk_pe(risk_count) = (((10 - R_o(ktrial, kcond) + 10 - o(ktrial,kcond)) - EV(ktrial, o(ktrial, kcond)+1)) ^2) - risk(risk_count) ;
                %                 end
                sspe                = (((10 - offers) + 10 - o(ktrial,kcond)) - EV(ktrial, o(ktrial, kcond)+1)).^2;
                risk(risk_count)    = sum(rPA .* sspe);
                
                risk_pe(risk_count) = (((10 - R_o(ktrial, kcond) + 10 - o(ktrial,kcond)) - EV(ktrial, o(ktrial, kcond)+1)) ^2) - risk(risk_count) ;
        end
        
        
        pc(ktrial, :)       = exp(beta1.*EV(ktrial,:)) ./ sum(exp(beta1.*EV(ktrial,:)));   % multinomial choice function
        lik(ktrial,kcond)   = pc(ktrial, o(ktrial,kcond)+1);
        EPc(ktrial)         = PA(ktrial,o(ktrial,kcond)+1);
        
        switch nmodel
            
            case 1
                PE(ktrial)       = r(ktrial, kcond) - EPc(ktrial);
                a_t(ktrial+1)   = a_t(ktrial) + 10.*lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 2
                PE(ktrial)       = r(ktrial,kcond) - EPc(ktrial);
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + 10.*lr2.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
            case 3
                PE(ktrial)       = (r(ktrial, kcond) - EPc(ktrial)).* reward;
                a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                
            case 4
                PE(ktrial)       = (r(ktrial, kcond) - EPc(ktrial)).*reward;
                if PE(ktrial)>0
                    a_t(ktrial+1)  = a_t(ktrial) + lr1.*PE(ktrial);
                else
                    a_t(ktrial+1)  = a_t(ktrial) + lr2.*PE(ktrial) ;                % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
                end
        end
    end
    all_pc(pc_counter:pc_counter+ntrial-1, :) = pc; % to store everything in a single array
    all_PE(pc_counter:pc_counter+ntrial-1, :) = PE; % to store everything in a single array
    all_EV(pc_counter:pc_counter+ntrial-1, :) = EV; % to store everything in a single array
    all_PA(pc_counter:pc_counter+ntrial-1, :) = PA; % to store everything in a single array
    
    pc_counter = pc_counter + ntrial;
end

negLL = - sum(log(lik(:)));

% lik = -lik;                                                                               % loglikelyhood vector
