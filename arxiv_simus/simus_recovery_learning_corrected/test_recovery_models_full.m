clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims  = 40;                           % nsubs to simulates
n_trial = 12;                           % ntrial per cond per session
n_sess  = 2;                            % nsession
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));% parameters of the simulation

modelspace = [1 2 3 4];
nfpm=[2 3 2 3];
nmods = numel(modelspace);

% set up conditions and mutliple sessions
%------------------------------------------
cond2learn  = -[12,9,6,3,0];
nc          = numel(cond2learn);
Ra          = repmat(cond2learn,1,n_sess);            % predator true accepance thereshold (logit intercept)
Rb          = repmat(3*ones(1,nc),1,n_sess);              % predator true accpetance moise (logit slope)
n_cond      = size(Ra,2);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 2.9+.2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)

Px_rnd          = .5+2.5*rand(n_sims,1);         %  Proposer  rating temperature
% Px_rnd          = 3+3*rand(n_sims,1);         %  Proposer  rating temperature
Plr1_rnd        = rand(n_sims,1);           %  Proposer  learning rate
Plr2_rnd        = rand(n_sims,1);           %  Proposer  learning rate


% setup estimation
%---------------------
options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);
parameters  = NaN(n_sims,3,nmods,nmods);        parametersLPP  = NaN(n_sims,3,nmods,nmods);
ll          = NaN(n_sims,nmods,nmods);        LPP            = NaN(n_sims,nmods,nmods);

% Sim loop
for ktm = modelspace  % ktm = k true model
    %----------
    for k_sim = 1:n_sims
        
        % pre-allocate
        O_mat = NaN(n_trial,n_cond);
        D_mat = NaN(n_trial,n_cond);
        
        % get params
        a0  = Pa_rnd(k_sim);
        b0  = Pb_rnd(k_sim);
        bX  = Px_rnd(k_sim);
        lr1 = Plr1_rnd(k_sim);
        lr2 = Plr2_rnd(k_sim);
        
        [O,D] = learning_models_timeseries([bX,lr1,lr2],[Ra;Rb],n_trial,a0,b0,ktm);
        
        
        lb = [0 0 0];          LB = [0 0 0];
        ub = [15 1 1];         UB = [Inf 10 1];
        ddb = ub - lb;
        
        for kem = modelspace
            
            n_rep           = 5;
            parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
            ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
            
            for k_rep = 1:n_rep
                x0 = lb + rand(1,3).*ddb;
                % %standard estimation
                [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                % %lalace estimation
                [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
            end
            [~,pos] = min(ll_rep);
            parameters(k_sim,:,ktm,kem)    =   parameters_rep(pos(1),:);
            ll(k_sim,ktm,kem)              =   ll_rep(pos(1),:);
            
            [~,posLPP] = min(LPP_rep);
            parametersLPP(k_sim,:,ktm,kem)      =   parametersLPP_rep(posLPP(1),:);
            LPP(k_sim,ktm,kem)                  =   LPP_rep(posLPP(1),:);
        end
    end
end



save('ML_recovery_2017_10_19_2')

