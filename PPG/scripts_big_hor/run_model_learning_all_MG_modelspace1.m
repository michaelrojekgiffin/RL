% this runs the first tier of the model comparison, comparing our 4 main
% models (2 choice prediction and 2 reward prediction, both with 1 and 2
% prediction errors).

clear
close all force
clc

%% find paths
cur_dir = pwd;
%  project_name = 'RL_PreyPredator';
project_name = 'RL/PPG'; % for use in michael's dropbox
study_name   = 'matlab';

findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,['data_',study_name]);

fl_list = dir(strcat(data_dir,filesep,'*DATA*.mat*'));
nsub = length(fl_list);

options     = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 1000000);

%% pre allocate
good_sub = NaN(nsub,1);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:10;
endow  = 10*ones(1,numel(offers));% parameters of the simulation

ntr = 60;
% pre-allocate
pred_offer_mat      = NaN(nsub, ntr);
pred_accept_mat     = NaN(nsub, ntr);
pred_reward_mat     = NaN(nsub, ntr);

prey_offer_mat      = NaN(nsub, ntr);
prey_accept_mat     = NaN(nsub, ntr);
prey_reward_mat     = NaN(nsub, ntr);
for k_sub = 1:nsub
    
    flnm = fullfile(data_dir,fl_list(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
            pred_offer_mat(k_sub, :)   =  data(:, 5)';  % column 5 is playerselect
            pred_accept_mat(k_sub, :)  =  data(:, 11)'; % column 11 is payoff
            pred_reward_mat(k_sub, :)  =  data(:, 9)';  % column 9 is payoff
            
            O                          = squeeze(pred_offer_mat(k_sub,:));
            D                          = squeeze(pred_accept_mat(k_sub,:));
            R_o                        = data(:, 6)';   % column 6 is offer of opponent
            
            role_dim                   = 1;             % tells which dimension to store loglik and parameters
        case 'prey'
            prey_offer_mat(k_sub, :)   =  data(:, 5)';  % column 5 is playerselect
            prey_accept_mat(k_sub, :)  =  data(:, 11)'; % column 11 is payoff
            prey_reward_mat(k_sub, :)  =  data(:, 9)';  % column 9 is payoff
            
            O                          = squeeze(prey_offer_mat(k_sub,:));
            D                          = squeeze(prey_accept_mat(k_sub,:));
            R_o                        = data(:, 6)';   % column 6 is offer of opponent
            
            role_dim                   = 2;             % tells which dimension to store loglik and parameters
    end
    
    
    for nmodel =1:4
        
        fprintf('Running model %d on subject %d of %d\n', nmodel, k_sub, nsub);
        
        n_rep           = 10;
        parameters_rep  = NaN(n_rep,5);     parametersLPP_rep  = NaN(n_rep,5);
        ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
        
        lb = [0 -15 0 0 0];          LB = [0 -Inf 0 0 0];
        ub = [5 0 5 1 1];         UB = [Inf Inf Inf 1 1];
        ddb = ub - lb;
        
        for k_rep = 1:n_rep
            x0 = lb + rand(1,5).*ddb;
            %standard estimation
            [parameters_rep(k_rep,1:5),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,nmodel, role, R_o),x0,[],[],[],[],lb,ub,[],options);
            %lalace estimation
            [parametersLPP_rep(k_rep,1:5),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,nmodel, role, R_o),x0,[],[],[],[],LB,UB,[],options);
        end
        
        [~,pos] = min(ll_rep);
        LEARN_parameters(k_sub,:,nmodel,role_dim)         =   parameters_rep(pos(1),:);
        LEARN_ll(k_sub,nmodel,role_dim)                   =   ll_rep(pos(1),:);
        
        [~,posLPP] = min(LPP_rep);
        LEARN_parametersLPP(k_sub,:,nmodel,role_dim)      =   parametersLPP_rep(posLPP(1),:);
        LEARN_LPP(k_sub,nmodel,role_dim)                  =   LPP_rep(posLPP(1),:);
    end
    
    
end


save('Learning_modelspace1_2017_09_03')