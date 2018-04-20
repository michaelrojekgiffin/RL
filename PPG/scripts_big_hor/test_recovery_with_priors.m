% this script is an adaptation of Mael's script that runs simulations for 4
% different models, choice prediction and reward prediction, and assesses
% recoverability and identifiability of each. It uses both Laplace and
% normal estimation. What I need to do is to use it in order to fit models
% to the PPG data with 2 learning rates (still need to work out exactly why
% that is), and with a risk parameter. This is going to produce lots of
% models (4 if I exclude models that only have 1 learning rate, which Mael
% recommends, choice, choice + risk, reward, reward + risk).
clear
close all
clc

project_name = '~/Dropbox/RL/PPG/'; % for use in michael's dropbox


% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
predprey_array  =  {'predator', 'prey'};
n_sims          = 30;                           % nsubs to simulates
n_trial         = 60;                           % ntrial per cond per session
n_sess          = 1;                            % nsession
offers          = 0:1:10;
suboffers       = 0:1:10;
endow           = 10*ones(1,numel(offers));     % parameters of the simulation
subendow        = 10*ones(1,numel(suboffers));  % parameters of the simulation
nmodel_array    = 1:4;                          % all the models to loop through
lr1_upper_bound  = 1;                            % this is the upper bound on the first learning rate, can be anywere between 0 and 11
lr2_upper_bound  = 1;                            % this is the upper bound on the first learning rate, can be anywere between 0 and 11

% set up conditions and mutliple sessions
%------------------------------------------
cond2learn  = -[5.3];%-[12,9,6,3,0];                        % all the intercepts the player needs to learn
nc          = numel(cond2learn);
Ra          = repmat(cond2learn,1,n_sess);          % predator true accepance thereshold (logit intercept)
% Rb          = repmat(3*ones(1,nc),1,n_sess);        % predator true accpetance noise (logit slope)
Rb          = repmat(.5*ones(1,nc),1,n_sess);        % predator true accpetance noise (logit slope)
n_cond      = size(Ra,2);
% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
% logitp = @(b,x) b(3)/2+(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));
%logitp = @(b,x)(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));
% Ubch = .10;

% Generate params
%-------------------
% priors I'm using the empirical ones estimated from subject data
% Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
% Pb_rnd          = 1.5+1*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
% Pb_rnd          = 2.5+1*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)

Px_rnd          = .5+1*rand(n_sims,1);% .5+2.5*rand(n_sims,1);         %  Proposer  rating temperature
% Px_rnd          = 3+3*rand(n_sims,1);         %  Proposer  rating temperature

Plr1_rnd        = lr1_upper_bound*rand(n_sims,1);      %  Proposer  learning rate
Plr2_rnd        = lr2_upper_bound*rand(n_sims,1);      %  Proposer  learning rate
% Plr2_rnd        = 1-Plr1_rnd;                         %  Proposer  learning rate - seeing if I force the two learning rates apart if I get better recovery

% PreAllocate
%---------------
A_mat = zeros(n_sims,n_trial,nc);
PE_mat = zeros(n_sims,n_trial,nc);
PAsub_mat = NaN(n_sims,numel(suboffers));

% setup estimation
%---------------------
options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000, 'Display', 'off');
parameters  = NaN(n_sims,3);        parametersLPP  = NaN(n_sims,3);
ll          = NaN(n_sims,1);        LPP            = NaN(n_sims,1);


% pre-allocate for the confusion matrix
%------------------------------------------
% confusion_parameters and confusion_parameresLPP are each mXnXpXz matrices
% have the following dimensions
%
% dim 1: predator and prey
% dim 2: simulated subject (so each entry is one simulated subject)
% dim 3: paraemters - so each entry is a different parameter of the model
% dim 4: the actual model we used to generate the simulated data 
% dim 5: the parameters estimated from 
con_parameters      = NaN(length(predprey_array), n_sims,3,4,4);        
con_parametersLPP   = NaN(length(predprey_array), n_sims,3,4,4);
con_ll              = NaN(length(predprey_array), n_sims,4,4);        
con_LPP             = NaN(length(predprey_array), n_sims,4,4);


% Sim loop
%----------
for predprey_count = 1:length(predprey_array)
    predprey = predprey_array{predprey_count};
    switch predprey
        case 'predator'
            priors   = load([project_name, 'data_priors/predator_priors.mat']);
        case 'prey'
            priors   = load([project_name, 'data_priors/prey_priors.mat']);
    end
    
    for model_iter = 1:length(nmodel_array)
        nmodel = nmodel_array(model_iter);
        
        for k_sim = 1:n_sims
            fprintf('running %s sim %d of %d, model %d of %d\n', predprey, k_sim, n_sims, model_iter, length(nmodel_array)); 
            % pre-allocate
            O_mat = NaN(n_trial,n_cond);
            D_mat = NaN(n_trial,n_cond);
            
            % get params
%             a0  = Pa_rnd(k_sim);
%             b0  = Pb_rnd(k_sim);
            a0  = priors.parameters(1);
            b0  = priors.parameters(2);
            
            bX  = Px_rnd(k_sim);
            lr1 = Plr1_rnd(k_sim);
            lr2 = Plr2_rnd(k_sim);
            
            % simulation
            %---------------------------------------------------------------------
            % O is offer
            % D is the reward, or decision outcome, 0 for loss, and 1 for win
            % PE is prediciton error
            % at is prediction of the intercept of the opponent's choice function
            %---------------------------------------------------------------------
            % for models 1 and 3, only 1 learning rate is actually used to
            % generate the data
            [O,D,PE,at, pd, R_o] = learning_models_timeseries_with_priors_MG([bX,lr1,lr2],[Ra;Rb],n_trial,a0,b0,nmodel, predprey);
            
            n_rep           = 5;
            
            % this loop allows me to fit all the different models to data
            % that is created from all different models, for confusion
            % matrix
            for con_model_iter = 1:length(nmodel_array)
                con_nmodel = nmodel_array(con_model_iter);

                if con_nmodel == 1 || con_nmodel == 3
                    numfreeparams = 2;
                    ub = [5 lr1_upper_bound];         UB = [Inf lr1_upper_bound];
                else
                    numfreeparams = 3;
                    ub = [5 lr1_upper_bound lr2_upper_bound];         UB = [Inf lr1_upper_bound lr2_upper_bound];
                end
                
                parameters_rep  = NaN(n_rep,numfreeparams);     parametersLPP_rep  = NaN(n_rep,numfreeparams);
                ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
                
                lb = zeros(1, numfreeparams);          LB = zeros(1, numfreeparams);
                
                ddb = ub - lb;
                
                for k_rep = 1:n_rep
                    x0 = lb + rand(1,numfreeparams).*ddb;
                    %standard estimation
                    [parameters_rep(k_rep,1:numfreeparams),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim_with_priors_MG(x,O,D,a0,b0,con_nmodel, predprey, R_o),x0,[],[],[],[],LB,UB,[],options);
                    % %                 [parameters_rep(k_rep,1:numfreeparams),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim_1lr_2017_09_26(x,O,D,a0,b0,nmodel, predprey),x0,[],[],[],[],LB,UB,[],options);
                    
                    %laplace estimation
                    [parametersLPP_rep(k_rep,1:numfreeparams),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning_with_priors_MG(x,O,D,a0,b0,con_nmodel, lr1_upper_bound, lr2_upper_bound, predprey, R_o),x0,[],[],[],[],LB,UB,[],options);
                end
                
                [~,pos]                              =   min(ll_rep);
                parameters(k_sim,1:numfreeparams)    =   parameters_rep(pos(1),:);
                ll(k_sim)                            =   ll_rep(pos(1),:);
                
                [~,posLPP]                           =   min(LPP_rep);
                parametersLPP(k_sim,1:numfreeparams) =   parametersLPP_rep(posLPP(1),:);
                LPP(k_sim)                           =   LPP_rep(posLPP(1),:);
                
                con_parameters(predprey_count, k_sim, 1:numfreeparams, nmodel, con_nmodel)        = parameters(k_sim,1:numfreeparams);
                con_parametersLPP(predprey_count, k_sim, 1:numfreeparams, nmodel, con_nmodel)     = parametersLPP(k_sim,1:numfreeparams);  
                con_ll(predprey_count, k_sim, nmodel, con_nmodel)                                 = ll(k_sim);
                con_LPP(predprey_count, k_sim, nmodel, con_nmodel)                                = LPP(k_sim);
                
                
            end
        end
        
    end
end


ws_nm = ['test_recovery_workspace_' num2str(n_sims), '_', date];
save(ws_nm)