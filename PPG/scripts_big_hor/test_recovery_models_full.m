clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims  = 60;                           % nsubs to simulates (should make it equal to number of subjects that I'm fitting to)
% importantly this must be divisible by 2 because I'm running this for both
% predator and prey
n_trial = 60;                           % ntrial per cond per session
n_sess  = 1;                            % nsession
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));% parameters of the simulation

role_array = ([ones(n_sims/2, 1); repmat(2, n_sims/2, 1)]);

modelspace = [1 2 3 4];
nfpm=[4 5 4 5];
nmods = numel(modelspace);

% set up conditions and mutliple sessions
%------------------------------------------
cond2learn  = -[6];
nc          = numel(cond2learn);
Ra          = repmat(cond2learn,1,n_sess);            % predator true accepance thereshold (logit intercept)
Rb          = repmat(3*ones(1,nc),1,n_sess);              % predator true accpetance moise (logit slope)
n_cond      = size(Ra,2);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
Pa_rnd          = -3 - 5*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = .4+2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
Px_rnd          = .05+3*rand(n_sims,1);         %  Proposer  rating temperature
Plr2_rnd        = 1*rand(n_sims,1);           %  Proposer  learning rate
Plr1_rnd        = 1*rand(n_sims,1);           %  Proposer  learning rate

% setup estimation
%---------------------
options     = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 1000000);

% all have an extra dimension added for predator and prey
parameters  = NaN(n_sims,5,nmods,nmods, 2);        parametersLPP  = NaN(n_sims,5,nmods,nmods, 2);
ll          = NaN(n_sims,nmods,nmods, 2);        LPP            = NaN(n_sims,nmods,nmods, 2);

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
        
        switch role_array(k_sim)
            case 1
                role = 'predator';
            case 2
                role = 'prey';
        end
        
        [O,D, opponent_o] = learning_models_timeseries([bX,a0,b0,lr1,lr2],[Ra;Rb],n_trial,ktm, role);
        
        lb = [0 -15 0 0 0];          LB = [0 -Inf 0 0 0];
        ub = [5 0 5 1 1];         UB = [Inf Inf Inf 1 1];
        ddb = ub - lb;
        
        for kem = modelspace
            
            fprintf('running %s true model %d, estimated model %d, for simulation %d of %d\n', role, ktm, kem, k_sim, n_sims);
            
            n_rep           = 10;
            parameters_rep  = NaN(n_rep,5);     parametersLPP_rep  = NaN(n_rep,5);
            ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
            
            for k_rep = 1:n_rep
                x0 = lb + rand(1,5).*ddb;
                % %standard estimation
                [parameters_rep(k_rep,1:5),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,kem, role, opponent_o),x0,[],[],[],[],LB,UB,[],options);
                % %laplace estimation
                [parametersLPP_rep(k_rep,1:5),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,kem, role, opponent_o),x0,[],[],[],[],LB,UB,[],options);
            end
            [~,pos] = min(ll_rep);
            parameters(k_sim,:,ktm,kem, role_array(k_sim))    =   parameters_rep(pos(1),:);
            ll(k_sim,ktm,kem, role_array(k_sim))              =   ll_rep(pos(1),:);
            
            [~,posLPP] = min(LPP_rep);
            parametersLPP(k_sim,:,ktm,kem, role_array(k_sim))      =   parametersLPP_rep(posLPP(1),:);
            LPP(k_sim,ktm,kem, role_array(k_sim))                  =   LPP_rep(posLPP(1),:);
        end
    end
end

save(['recovery', num2str(n_sims) 'sims_2018_12_03'])

