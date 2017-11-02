% this script fits the model to subject data, and then runs a simulation
% for each subject using the parameters estimated from their data in order
% to compare
% clear
% close all
% clc

cd ~/Dropbox/RL/PPG/
%=========================================================================
%-------------------------------------------------------------------------
% EDIT - this part needs to be tweaked, depending on what I'm looking for
%-------------------------------------------------------------------------

% nsub          = 166;   % hor dataset - number of plots to make, maximum is 166 for hor dataset (from data_matlab) (i.e. length(fl_dir) )
% nsub          = 48;   % old_hor dataset - number of plots to make, maximum is 48 (i.e. length(fl_dir) )
% nsub          = 50;   % ital dataset - number of plots to make, maximum is 50 (i.e. length(fl_dir) )
nsub          = 10;   % OT data - number of subjects is 27, and each playes in both roles for a total of 54
plot_ind      = false;
nsim          = 15;
plot_all_data = false;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%=========================================================================

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% changes for each dataset
% Importantly, at the end of the first cell there's code that saves all the
% important variables to the specific data directories of each experiment,
% that needs to be changed
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
cur_dir      = pwd;
% data_dir     = fullfile(cur_dir,'data_matlab');             % hor dataset
% data_dir     = fullfile(cur_dir,'data_old_hor');             % old_hor dataset
% data_dir     = fullfile(cur_dir,'data_ital');             % ital dataset
data_dir     = fullfile(cur_dir,'data_OT');             % OT dataset
fl_dir       = dir(strcat(data_dir,filesep,'DATA_sub*'));

% hor and old_hor datasets have the same columns for these
% sub_o_col       = 5;
% sub_r_col       = 11;
% opponent_o_col  = 6;

% % for ital dataset, these are the columns
% sub_o_col       = 6;
% sub_r_col       = 12;
% opponent_o_col  = 7;

% for OT dataset, these are the columns. Importantly, in this dataset
% subjects have multiple "sessions", meaning they learn the same
% distribution over multiple blocks that each must begin at 0
sub_o_col       = 7;
sub_r_col       = 96;
nc              = 1;  % number of different "conditions", i.e. separate distributions subjects played against (just 1)
n_trial_per_cond= 20; % NOT TO BE CONFUSED WITH ntrial, WHICH IS TOTAL NUMBER OF TRIALS, THIS IS NUMBER OF TRIALS PER CONDITION
n_sess           = 3; % number of different blocks subjects played within each role
% opponent_o_col  = 7; % changes depending on the role, is specified near
% begininning of first loop for OT dataset

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

% define function
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
nmodel_array    = 1:4;                          % all the models to loop through
lr_upper_bound  = 11;                            % this is the upper bound on the first learning rate, can be anywere between 0 and 11

options      = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000);

parametersLPP  = NaN(nsim,3);
LPP            = NaN(nsim,1);

k_prey       = 0;
k_pred       = 0;

pred_ameters        = NaN(nsub, 3, length(nmodel_array)); % parameters of predator model for each subject playing as predator
prey_ameters        = NaN(nsub, 3, length(nmodel_array)); % parameters of prey model for each subject playing as prey
prey_freq           = []; % prey frequencies
pred_freq           = []; % predator frequencies
pred_dist           = NaN(nsub, 60, length(nmodel_array)); % predator distribution, 60 because 60 trials
prey_dist           = NaN(nsub, 60, length(nmodel_array)); % prey distribution
pred_sim_dist       = NaN(nsub, 60); % distribution created from the simulations of predator (for each subject)
prey_sim_dist       = NaN(nsub, 60); % distribution created from the simulations of prey (for each subject)
pred_opponent_probs = []; % distribution each predator played against
prey_opponent_dist  = []; % distribution each prey played against

                
EV_sub              = NaN(nsub, 60, length(nmodel_array));             % expected value of each decision given parameters recovered from data (each subject has own row, each trial has own column)
PA_sub              = NaN(nsub, length(offers), length(nmodel_array)); % probaility that each investment will succeed given params of model
V_sub               = NaN(nsub, length(nmodel_array));              % posterior intercept of subject decision function given params of model
EV_sub_posterior    = NaN(nsub, length(offers), length(nmodel_array)); % expected value of every possible offer with posterior intercept given recovered model params
pc_sub              = NaN(nsub, 60, length(offers), length(nmodel_array));

BIC                 = NaN(nsub, length(nmodel_array)); %
LPP                 = NaN(nsub, length(nmodel_array)); %
prey_BIC            = NaN(nsub, length(nmodel_array)); %
prey_LPP            = NaN(nsub, length(nmodel_array)); %
pred_BIC            = NaN(nsub, length(nmodel_array)); %
pred_LPP            = NaN(nsub, length(nmodel_array)); %

PE_sub              = NaN(nsub, 60, length(nmodel_array));             % prediction error for each subject on each trial for each model
pred_PE             = NaN(nsub, 60, length(nmodel_array));             % prediction error for each subject on each trial for each model
prey_PE             = NaN(nsub, 60, length(nmodel_array));             % prediction error for each subject on each trial for each model

pred_EV             = NaN(nsub, 60, length(nmodel_array));             % EV
prey_EV             = NaN(nsub, 60, length(nmodel_array));             

risk_sub              = NaN(nsub, 60, length(nmodel_array));             % risk prediction for each subject and each model
risk_sub_pe           = NaN(nsub, 60, length(nmodel_array));             % risk prediction error for each subject on each trial for each model

pred_risk             = NaN(nsub, 60, length(nmodel_array));            
pred_risk_pe          = NaN(nsub, 60, length(nmodel_array)); 

prey_risk             = NaN(nsub, 60, length(nmodel_array));            
prey_risk_pe          = NaN(nsub, 60, length(nmodel_array)); 
for k_sub    = 1:nsub
    
    flnm     = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
            k_pred   = k_pred+1;
            predprey = 'predator';
            opponent = 'prey';
            priors   = load('data_priors/predator_priors.mat');
            opponent_o_col  = 9; % changes depending on the role
        case 'prey'
            k_prey   = k_prey+1;
            predprey = 'prey';
            opponent = 'predator';
            priors   = load('data_priors/prey_priors.mat');
            opponent_o_col  = 8; % changes depending on the role
    end
    %
    for k_model = 1:length(nmodel_array)
        % have to keep track of the models
        
        nmodel = nmodel_array(k_model);
        
        %         fprintf('estimating for subject %s, model %d, %s, %d out of %d\n', fl_dir(k_sub).name(6:11), nmodel, predprey, k_sub, nsub);
        fprintf('estimating for subject %s, model %d, %s, %d out of %d\n', fl_dir(k_sub).name(10:15), nmodel, predprey, k_sub, nsub);
        
        % this is true for each dataset in which subjects played either in
        % only 1 role or in 1 role per session (so the 3 hormone datasets)
%         sub_o            = data(:,sub_o_col);       % subject offer
%         sub_r            = data(:,sub_r_col);       % subject win/lose (logical)
%         opponent_o       = data(:,opponent_o_col);  %  choice of the opponent
        
        % this only applies for datasets with multiple sessions, in the
        % case of OT, it has 3 for both predator and prey.
        ses_counter = 1;
        for SES = 1:3 % because OT has 3 sessions per role
            
            sub_o(1:20,SES)            = data(ses_counter:ses_counter+19,sub_o_col);       % subject offer
            sub_r(1:20,SES)            = data(ses_counter:ses_counter+19,sub_r_col);       % subject win/lose (logical)  
            opponent_o(1:20,SES)       = data(ses_counter:ses_counter+19,opponent_o_col);  %  choice of the opponent
            ses_counter                = ses_counter+20;
        end
        
        
        %%%%%%%% for some reason this freaked out when I did it with the OT
        %%%%%%%% data
%         ntrial           = length(data);
        ntrial           = size(data, 1);
        
        n_rep            = 10;
        parameters_rep   = NaN(n_rep,2);
        ll_rep           = NaN(n_rep,1);
        
        
        if nmodel == 1 || nmodel == 3
            numfreeparams = 2;
            ub = [5 lr_upper_bound];         UB = [Inf lr_upper_bound];
        else
            numfreeparams = 3;
            ub = [5 lr_upper_bound 1];         UB = [Inf lr_upper_bound 1];
        end
        
        parametersLPP_rep  = NaN(n_rep,numfreeparams);
        LPP_rep          = NaN(n_rep,1);
        
        lb = zeros(1, numfreeparams);          LB = zeros(1, numfreeparams);
        
        ddb = ub - lb;
        
        for k_rep = 1:n_rep
            x0 = lb + rand(1,numfreeparams).*ddb;
            %laplace estimation
            [parametersLPP_rep(k_rep,1:numfreeparams),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2_MG(x,sub_o,sub_r,priors.parameters(1),priors.parameters(2),nmodel, lr_upper_bound, predprey, opponent_o),x0,[],[],[],[],LB,UB,[],options);
        end
% %         
% %         params = [2 .3]
% %         o = sub_o
% %         r = sub_r
% %         ,r,a0,b0,nmodel, predprey, R_o
        
        
        [~,posLPP]                            =   min(LPP_rep);
        parametersLPP(k_sub,1:numfreeparams)  =   parametersLPP_rep(posLPP(1),:);
        
        LPP(k_sub, k_model)                   =   LPP_rep(posLPP(1),:);
        
        [~, BIC(k_sub, k_model)]              =   aicbic(-LPP(k_sub, k_model), numfreeparams, ntrial);
        
        
        % record subject's expected values (sub_EV), estimates of probability
        % of success (sub_PA), and intercept (sub_V), given the parameters
        % estimated above
        %     [~, all_EV, all_PA, all_V] = ModelEstimation_2params_2017_09_11(parameters(k_sub,:), priors.parameters(1:2),sub_o,sub_r,1,predprey);
        
        [~, all_EV, all_PA, all_V, all_pc, all_PE, all_risk, all_risk_pe]  = laplace_priors_learning2_MG(parametersLPP(k_sub,1:numfreeparams),sub_o,sub_r,priors.parameters(1),priors.parameters(2), nmodel, lr_upper_bound, predprey, opponent_o);
        V_sub(k_sub, k_model)               = all_V(end);     % posterior intercept
        PA_sub(k_sub, :, k_model)           = all_PA(end, :); % posterior prob of success
        EV_sub_posterior(k_sub, :, k_model) = all_EV(end, :); % posterior EV of all offers
        pc_sub(k_sub, :, :, k_model)        = all_pc;         % expected offer, i.e. probability of each offer being selected on each trial
        PE_sub(k_sub, :, k_model)           = all_PE;         % prediction error
        
        risk_sub(k_sub, :, k_model)         = all_risk;       % expected risk
        risk_sub_pe(k_sub, :, k_model)      = all_risk_pe;    % risk prediction error
        for ee = 1:length(all_EV)
            EV_sub(k_sub, ee, k_model) = all_EV(ee, sub_o(ee)+1); % EV of each offer selected by subject
        end
        
        
        % Get slope and intercept of distribution against which the subject is
        % playind (i.e. estimate their true "acceptance function")
        opponent_o_freq = zeros(1, 11);
        succes_probs = zeros(length(opponent_o_freq), 1); % probability of success
        for k = 0:10
            opponent_o_freq(k+1)       = (sum(opponent_o(:)==k)) / length(opponent_o);
            switch predprey
                case 'prey'
                    succes_probs(k+1) = sum(opponent_o_freq(1:k+1));
                    opponent_role = 'predator';
                case 'predator'                                         % predator only wins if prey invests less
                    succes_probs(k+1) = sum(opponent_o_freq(1:k));     % sum all probabilities below chosen offer
                    opponent_role = 'prey';
            end
        end
        empirical_EV = zeros(length(offers), 1);
        switch predprey
            case 'prey'
                empirical_EV = (endow - offers).* succes_probs';
            case 'predator'
                empirical_EV = (endow - offers) + ((endow - offers).* succes_probs');
        end
        
        % here I'm fitting to the opponent data in order to get the logit
        % choice function that subjects were playing against in order to
        % use this in the simulation below.
        n_rep            = 10;
        parameters_rep   = NaN(n_rep,3);
        ll_rep           = NaN(n_rep,1);
        opponent_o_array = reshape(opponent_o, [], 1);
        
        for k_rep = 1:n_rep
            [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) laplace_priors_priors_MG_2017_10_03(x,opponent_o_array, opponent_role),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
        end
        [~,pos]             = min(ll_rep);
        opponent_parameters = parameters_rep(pos(1),:);
        ll                  = ll_rep(pos(1),:);
        
        
        
        % initialize for simulation
        B0      = parametersLPP(k_sub, 1);   % proposer temperature, estimated from subject data (specific to each subject)
        a0      = parametersLPP(k_sub, 2);   % proposer learning rate, estimated from subject data (specific to each subject)
        v0      = priors.parameters(1);      % proposer initial prior on threshold (estimated from concatenated first trial investments)
        B       = priors.parameters(2);      % proposer noise (proposer's estimate of slope of opponent, (estimated from concatenated first trial investments))
        
        % pre-allocat,
        %--------------------------
        PA      = NaN(ntrial,numel(offers), nsim); % estimated probability of accepting all offers
        EV      = NaN(ntrial,numel(offers), nsim); % estimated expected value of all offers
        V       = NaN(ntrial+1, nsim);             % estimated opponent intercept
        kO      = NaN(ntrial, nsim);               % selected offer, in the 1:numel(offer) space
        O       = NaN(ntrial, nsim);               % selected offer, in the euro offer spavce
        Pd      = NaN(ntrial, nsim);               % trial  probability of accepting the selected offer
        D       = NaN(ntrial, nsim);               % final decision of accepting the selected offer
        CPE     = NaN(ntrial, nsim);               % Choice Prediction error
        pc_sim  = NaN(ntrial, numel(offers), nsim);% expected offer
        % initialize simulation
        %--------------------------
        for i_sim = 1:nsim
            V(1, i_sim) = v0;
            
            % simulation
            %---------------------------------------------------------------------
            % input
            % bX is the temperature
            % lr1 is learning rate 1
            % lr2 is learing rate 2 (only actually used for models 2 and 4)
            % opponent_parameters were fitted from the distribution against
            % which subjects were playing, in the actual function they're
            % written as Ra and Rb, the intercept and slope of the
            % opponent's logic function, respectively
            % Ra is the true intercept of the opponent's logit function
            % Rb is the true slope of the opponent's logic functiuon
            % ntrial is the number of trials to simulate
            % v0 is the prior estimate of the intercept
            % B is the prior estimate of the slope
            % opponent_o is the distribution of the opponent's decisions
            %
            % output
            % O is offer
            % D is the reward, or decision outcome, 0 for loss, and 1 for win
            % PE is prediciton error
            % V is prediction of the intercept of the opponent's choice function
            % pd is posterior density distribution
            % R_o is the simlation of the opponent's behavior, given the
            % model and the known logit function used as input
            %
            %---------------------------------------------------------------------
            % for models 1 and 3, only 1 learning rate is actually used to
            % generate the data
            bX      =    parametersLPP(k_sub, 1);   % proposer temperature, estimated from subject data (specific to each subject)
            lr1     =    parametersLPP(k_sub, 2);   % proposer learning rate 1, estimated from subject data (specific to each subject)
            lr2     =    parametersLPP(k_sub, 3);   % proposer learning rate 2, estimated from subject data (specific to each subject, will be NaN for models 1 and 3)
            
            [O(:, i_sim),D(:, i_sim), CPE(:, i_sim),V(:, i_sim), pc_sim(:, :, i_sim), R_o, PA(:, :, i_sim), EV(:, :, i_sim)] = learning_models_timeseries_MG_2017_10_03([bX,lr1,lr2],opponent_parameters,size(sub_o, 1),v0,B,nmodel, predprey, opponent_o);
            
        end
        % average across all simulations
        V_sim              = mean(V(end, :)); % average of V (opponent's intercept) estimate across simulations
        % pre-allocate
        EV_sim_prior       = NaN(1,length(offers));
        EV_sim_posterior   = NaN(1,length(offers));
        PA_sim_prior       = NaN(1,length(offers));
        PA_sim_posterior   = NaN(1,length(offers));
        pc_sim_mean        = NaN(ntrial, length(offers));
        O_mean             = NaN(ntrial, 1);
        
        for m_sim = 1:length(offers)
            EV_sim_prior(m_sim)     = mean(EV(1, m_sim, :));      % row 1 is prior
            EV_sim_posterior(m_sim) = mean(EV(ntrial, m_sim, :)); % last row is posterior
            PA_sim_prior(m_sim)     = mean(PA(1, m_sim, :));
            PA_sim_posterior(m_sim) = mean(PA(ntrial, m_sim, :));
        end
        % take average of offers made by simulation, averages within trial,
        % across simulations
        for m_trial = 1:ntrial
            O_mean(m_trial)      = mean(O(m_trial, :));
            for m_offer = 1:numel(offers)
                pc_sim_mean(m_trial, m_offer) = mean(pc_sim(m_trial, m_offer, :));
            end
        end
        sub_o_array = reshape(sub_o, [], 1);
        % record all parameters and all subject data for posterity
        switch predprey
            case 'predator'
                pred_ameters(k_pred, :, k_model)        = parametersLPP(k_sub, :); % column 1 is temperature, column 2 and 3 learning rate
                pred_dist(k_pred, :, k_model)           = sub_o_array;
                pred_sim_dist(k_pred, :)                = O_mean;
                pred_opponent_probs(k_pred, :)          = succes_probs;
                
                pred_PE(k_pred, :, k_model)              = PE_sub(k_sub, :, k_model);   % prediction error for each subject on each trial for each model
                pred_EV(k_pred, :, k_model)              = EV_sub(k_sub, :, k_model);   % prediction error for each subject on each trial for each model
                
                pred_risk(k_pred, :, k_model)           = risk_sub(k_sub, :, k_model);   % risk prediction 
                pred_risk_pe(k_pred, :, k_model)        = risk_sub_pe(k_sub, :, k_model);   % risk prediction 
                
                pred_LPP(k_pred, k_model)               = LPP(k_sub, k_model);
                [~, pred_BIC(k_pred, k_model)]          = aicbic(-LPP(k_sub, k_model), numfreeparams, ntrial);
                
            case 'prey'
                prey_ameters(k_prey, :, k_model)        = parametersLPP(k_sub, :); % column 1 is temperature, column 2 and 3 learning rate
                prey_dist(k_prey, :, k_model)           = sub_o_array;
                prey_sim_dist(k_prey, :)                = O_mean;
                prey_opponent_probs(k_prey, :)          = succes_probs;
                
                prey_PE(k_prey, :, k_model)              = PE_sub(k_sub, :, k_model);   % prediction error for each subject on each trial for each model
                prey_EV(k_prey, :, k_model)              = EV_sub(k_sub, :, k_model);   % prediction error for each subject on each trial for each model
                
                
                pred_risk(k_prey, :, k_model)           = risk_sub(k_sub, :, k_model);   % risk prediction 
                pred_risk_pe(k_prey, :, k_model)        = risk_sub_pe(k_sub, :, k_model);   % risk prediction error
                
                
                prey_LPP(k_prey, k_model)               = LPP(k_sub, k_model);
                [~, prey_BIC(k_prey, k_model)]          = aicbic(-LPP(k_sub, k_model), numfreeparams, ntrial);
        end
 
        % plot everything, per individual if plot_ind = true
        if plot_ind
            if ~plot_all_data
                figure
                subplot(6, 3, [2:3, 5:6, 8:9])
                hold on
                plot(offers,succes_probs,'-r')
                plot(offers,logitp(priors.parameters(1:2), offers),'-.k', 'Linewidth', 1.5)
                plot(offers,logitp([V_sub(k_sub, k_model), B(end)],offers),'--b')
                plot(offers,logitp([V_sim, B],offers),'--g')
                legend('true \phi','prior \phi','posterior \phi', 'simulated posterior \phi')
                xlabel('Offer (euro)')
                ylabel('p(Accept)')
                title(['Subject ' num2str(k_sub) ' (' predprey ')' ' model ', num2str(k_model)]);
                
                % this plot is wrong, it is depicting simulation, should be
                % depciting subject - need to simply add a line for the subject
                subplot(6, 3, [11:12, 14:15, 17:18])
                hold on
                plot(offers, empirical_EV, '-r')
                plot(offers, EV_sim_prior, '-.k', 'Linewidth', 1.5)
                plot(offers, EV_sub_posterior(k_sub, :, k_model), '--b')
                plot(offers, EV_sim_posterior, '--g')
                legend('true', 'prior', 'posterior', 'simulated posterior')
                xlabel('Offer (euro)')
                ylabel('Expected Value')
                
                inv_counter = 1;
                fibb = 0;
                for ii = 1:ntrial/10
                    subplot(6, 3, ii + fibb)
                    fibb = fibb+2;
                    inv_row = str2double(sprintf('%d0', ii));
                    
                    tmp_pc_sub = squeeze(pc_sub(k_sub, inv_counter:inv_row,:, k_model));
                    
                    bar(0:10, tmp_pc_sub(end,:), 'b')
                    hold on
                    bar(-.5:1:9.5, pc_sim_mean(end,:), 'g')
                    plot(sub_o(inv_row), .1, 'Xr', 'MarkerSize', 10)
                    
                    if ii == 1
                        title('Real vs. simulated subject expected offer, every 10 trials')
                        legend('Subject', 'Simulation', 'subject selects')
                    end
                    hold off
                    inv_counter = inv_counter+10;
                end
                
            else
                figure
                subplot(6, 3, 1:3:7)
                plot(1:length(data), sub_o, 'b-o', ...
                    1:length(data), opponent_o, 'r--*',...
                    1:length(data), O_mean, '-.ks')
                legend(['Subject/' predprey], ['Opponent/' opponent], 'Simulated')
                title(['Subject (' predprey ') opponent (' opponent, ') and simiulation']);
                
                subplot(6, 3, 10:3:18)
                histogram(sub_o, 11)
                hold on
                histogram(O_mean, 11)
                title('Real vs. simulated subject historgram, all trials')
                legend('Subject', 'Simulation')
                
                inv_counter = 1;
                fibb = 0;
                for ii = 1:ntrial/10
                    subplot(6, 3, ii+ii + fibb)
                    fibb = fibb+1;
                    inv_row = str2double(sprintf('%d0', ii));
                    histogram(sub_o(inv_counter:inv_row,:), 11)
                    hold on
                    histogram(O_mean(inv_counter:inv_row,:), 11)
                    if ii == 1
                        title('Real vs. simulated subject historgram, every 10 trials')
                        legend('Subject', 'Simulation')
                    end
                    hold off
                    inv_counter = inv_counter+10;
                end
                
                subplot(6, 3, 3:3:9)
                hold on
                plot(offers,succes_probs','-r')
                plot(offers,logitp(priors.parameters(1:2), offers),'-.k', 'Linewidth', 1.5)
                plot(offers,logitp([V_sub(k_sub, k_model), B(end)],offers),'--b')
                plot(offers,logitp([V_sim, B],offers),'--g')
                legend('true \phi','prior \phi','posterior \phi', 'simulated posterior \phi')
                xlabel('Offer (euro)')
                ylabel('p(Accept)')
                title(['Subject (' predprey ')']);
                
                subplot(6, 3, 12:3:18)
                hold on
                plot(offers, empirical_EV, '-r')
                plot(offers, EV_sim_prior, '-.k', 'Linewidth', 1.5)
                % Here i have to divise a way to keep track of the model
                % that I'm on
                plot(offers, EV_sub_posterior(k_sub, :, k_model), '--b')
                plot(offers, EV_sim_posterior, '--g')
                legend('true', 'prior', 'posterior', 'simulated posterior')
                xlabel('Offer (euro)')
                ylabel('Expected Value')
                
            end
        end
    end
    
end


figure
imagesc(BIC)
colsize = size(BIC, 2);
colcounter = 0;
rowcounter = 1;
for ii = 1:numel(BIC)
    colcounter = colcounter+1;
    text(colcounter, rowcounter, num2str(BIC(rowcounter, colcounter)), 'Color', 'r')
    if colcounter == colsize
        colcounter = 0;
        rowcounter = rowcounter + 1;
    end
end
colormap default
title(sprintf('BIC'))
colorbar

% getting rid of NaN values in the role specifc BIC variables
pred_BIC(~any(~isnan(pred_BIC), 2), :) = [];
prey_BIC(~any(~isnan(prey_BIC), 2), :) = [];

% these are for hor
% save([pwd '/data_matlab/pred_ameters.mat'], 'pred_ameters');
% save([pwd '/data_matlab/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_matlab/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_matlab/pred_BIC.mat'], 'pred_BIC');
% save([pwd '/data_matlab/prey_BIC.mat'], 'prey_BIC');
% save([pwd '/data_matlab/pred_LPP.mat'], 'pred_LPP');
% save([pwd '/data_matlab/prey_LPP.mat'], 'prey_LPP');

% for old_hor
% save([pwd '/data_old_hor/pred_ameters.mat'], 'pred_ameters');
% save([pwd '/data_old_hor/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_old_hor/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_old_hor/pred_BIC.mat'], 'pred_BIC');
% save([pwd '/data_old_hor/prey_BIC.mat'], 'prey_BIC');
% save([pwd '/data_old_hor/pred_LPP.mat'], 'pred_LPP');
% save([pwd '/data_old_hor/prey_LPP.mat'], 'prey_LPP');

% for ital
% save([pwd '/data_ital/pred_ameters.mat'], 'pred_ameters');
% save([pwd '/data_ital/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_ital/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_ital/pred_BIC.mat'], 'pred_BIC');
% save([pwd '/data_ital/prey_BIC.mat'], 'prey_BIC');
% save([pwd '/data_ital/pred_LPP.mat'], 'pred_LPP');
% save([pwd '/data_ital/prey_LPP.mat'], 'prey_LPP');

% for OT
% save([pwd '/data_OT/pred_ameters.mat'], 'pred_ameters');
% save([pwd '/data_OT/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_OT/prey_ameters.mat'], 'prey_ameters');
% save([pwd '/data_OT/pred_BIC.mat'], 'pred_BIC');
% save([pwd '/data_OT/prey_BIC.mat'], 'prey_BIC');
% save([pwd '/data_OT/pred_LPP.mat'], 'pred_LPP');
% save([pwd '/data_OT/prey_LPP.mat'], 'prey_LPP');
% 
% save([pwd '/data_OT/PE_sub.mat'], 'PE_sub');
% save([pwd '/data_OT/EV_sub.mat'], 'EV_sub');
% save([pwd '/data_OT/prey_EV.mat'], 'prey_EV');
% save([pwd '/data_OT/prey_PE.mat'], 'prey_PE');
% save([pwd '/data_OT/pred_PE.mat'], 'pred_PE');
% save([pwd '/data_OT/pred_EV.mat'], 'pred_EV');
pred_LPP(~any(~isnan(pred_LPP), 2), :) = [];
prey_LPP(~any(~isnan(prey_LPP), 2), :) = [];


nfpm=[2 3 2 3]; % number of free parameters
bic = NaN(k_pred, length(nmodel_array));
for k_true = 1:length(nmodel_array)
    
    %     MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
    LL = pred_LPP(:,k_true);
    bic(:,k_true)=-2*-LL + nfpm(k_true)*log(nc*n_trial_per_cond*n_sess); % l2 is already positive
   
%     for k_est= 1:length(nmodel_array)
%         bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial_per_cond*n_sess); % l2 is already positive
%     end
%     [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
%     % [postBMC,outBMC]=VBA_groupBMC(-LL');
%     BMC_output(k_true).post = postBMC;
%     BMC_output(k_true).out = outBMC;
%     
%     Ep(k_true,:) = 100*BMC_output(k_true).out.ep;
end
[pred_postBMC,pred_outBMC]=VBA_groupBMC(-bic'./2);


nfpm=[2 3 2 3]; % number of free parameters
bic = NaN(k_prey, length(nmodel_array));
for k_true = 1:length(nmodel_array)
    
    %     MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
    LL = prey_LPP(:,k_true);
    bic(:,k_true)=-2*-LL + nfpm(k_true)*log(nc*n_trial_per_cond*n_sess); % l2 is already positive
   
%     for k_est= 1:length(nmodel_array)
%         bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial_per_cond*n_sess); % l2 is already positive
%     end
%     [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
%     % [postBMC,outBMC]=VBA_groupBMC(-LL');
%     BMC_output(k_true).post = postBMC;
%     BMC_output(k_true).out = outBMC;
%     
%     Ep(k_true,:) = 100*BMC_output(k_true).out.ep;
end
[prey_postBMC,prey_outBMC]=VBA_groupBMC(-bic'./2);


%% write to file
% write to file for analysis in R
% [H P] = ttest2(pred_ameters(2, :), prey_ameters(2, :))
% [H P] = ttest2(pred_ameters(1, :), prey_ameters(1, :))

% here I'm going to write everything, all the parameters, their liklihoods,
% and their BIC's, to a file for analysis in R. The numbers preceded by m
% next to the parameters and the liklihoods represent the model they
% correspond to. 
fitted_header = {'sub_name', 'role', 'model', 'alpha1', 'alpha2' 'beta', 'LL', 'BIC'};
fitted_length = nsub;
fitted_cell = cell(fitted_length*length(nmodel_array), length(fitted_header)); % each subject has one row for each model

sub_counter = 1;
sub_tc = 0; % subject trial counter
k_pred = 0;
k_prey = 0;
for ii = 1:nsub
    %     sub_tc   = sub_tc+1;
    flnm     = fullfile(data_dir,fl_dir(ii).name);
    load(flnm)
    %     % getting subject's name for hor
    %     if length(flnm) == 61
    %         sub_name = flnm(end-4);
    %     elseif length(flnm) == 62
    %         sub_name = flnm(end-5:end-4);
    %     elseif length(flnm) == 63
    %         sub_name = flnm(end-6:end-4);
    %     end
    
    % getting subject's name for old_hor
    %     switch role
    %         case 'predator'
    %             if length(flnm) == 72
    %                 sub_name = flnm(end-14:end-13);
    %             elseif length(flnm) == 71
    %                 sub_name = flnm(end-13);
    %             end
    %         case 'prey'
    %             if length(flnm) == 68
    %                 sub_name = flnm(end-10:end-9);
    %             elseif length(flnm) == 67
    %                 sub_name = flnm(end-9);
    %             end
    %     end
    
    %     % getting subject's name for ital
    %     if length(flnm) == 59
    %         sub_name = flnm(end-4);
    %     elseif length(flnm) == 60
    %         sub_name = flnm(end-5:end-4);
    %     end
    
    % getting subject's name for OT
    switch role
        case 'predator'
            if length(flnm) == 72
                sub_name = flnm(end-18:end-13);
%             elseif length(flnm) == 71
%                 sub_name = flnm(end-13);
            end
        case 'prey'
            if length(flnm) == 68
                sub_name = flnm(end-14:end-9);
%             elseif length(flnm) == 67
%                 sub_name = flnm(end-9);
            end
    end
    
    switch role
        case 'predator'
            k_pred   = k_pred+1;
            
            for mm = 1:length(nmodel_array)
                fitted_cell{sub_counter+sub_tc, 1} = sub_name;
                fitted_cell{sub_counter+sub_tc, 2} = 'predator';
                fitted_cell{sub_counter+sub_tc, 3} = mm;
                
                fitted_cell{sub_counter+sub_tc, 4} = pred_ameters(k_pred, 2, mm); % column 2 of pred_ameters is alpha (first learning rate)
                fitted_cell{sub_counter+sub_tc, 5} = pred_ameters(k_pred, 3, mm); % column 3 of pred_ameters is alpha (second learning rate)
                fitted_cell{sub_counter+sub_tc, 6} = pred_ameters(k_pred, 1, mm); % column 1 of pred_ameters is beta
                
                fitted_cell{sub_counter+sub_tc, 7} = pred_LPP(k_pred, mm);
                fitted_cell{sub_counter+sub_tc, 8} = pred_BIC(k_pred, mm);
                
                sub_tc                             = sub_tc + 1;
            end
            
        case 'prey'
            k_prey   = k_prey+1;
            for mm = 1:length(nmodel_array)
                
                fitted_cell{sub_counter+sub_tc, 1} = sub_name;
                fitted_cell{sub_counter+sub_tc, 2} = 'prey';
                fitted_cell{sub_counter+sub_tc, 3} = mm;
                
                fitted_cell{sub_counter+sub_tc, 4} = prey_ameters(k_prey, 2, mm); % column 2 of prey_ameters is alpha (first learning rate)
                fitted_cell{sub_counter+sub_tc, 5} = prey_ameters(k_prey, 3, mm); % column 3 of prey_ameters is alpha (second learning rate)
                fitted_cell{sub_counter+sub_tc, 6} = prey_ameters(k_prey, 1, mm); % column 1 of prey_ameters is beta
                
                fitted_cell{sub_counter+sub_tc, 7} = prey_LPP(k_prey, mm);
                fitted_cell{sub_counter+sub_tc, 8} = prey_BIC(k_prey, mm);
                
                sub_tc                             = sub_tc + 1;
            end
    end
    sub_tc      = 0;
    sub_counter = sub_counter + length(nmodel_array);
end
%
% %
% %
% fid = fopen('/Users/michaelgiffin/Carsten PhD/hormones/data/modeling/OT_fitted.txt', 'w');
% for ii = 1:length(fitted_cell)+1 % +1 for header
%     if ii == 1
%         fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', fitted_header{1,:});
%     else
%         fprintf(fid, '%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\n', fitted_cell{ii-1, :});
%     end
% end
% 
% fid = fopen('/Users/michaelgiffin/Carsten PhD/OT_data/RL/data/OT_fitted.txt', 'w');
% for ii = 1:length(fitted_cell)+1 % +1 for header
%     if ii == 1
%         fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', fitted_header{1,:});
%     else
%         fprintf(fid, '%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\n', fitted_cell{ii-1, :});
%     end
% end

%% Get all latent variables in table format
fitted_header = {'sub_name', 'role', 'model', 'EV', 'PE', 'risk', 'risk_pe'};
fitted_length = nsub*ntrial*nmodel*length(nmodel_array);
fitted_cell = cell(fitted_length, length(fitted_header)); % each subject has one row for each model

sub_counter = 1;
sub_tc = 0; % subject trial counter
k_pred = 0;
k_prey = 0;
for ii = 1:nsub
    flnm     = fullfile(data_dir,fl_dir(ii).name);
    load(flnm)
    %     % getting subject's name for hor
    %     if length(flnm) == 61
    %         sub_name = flnm(end-4);
    %     elseif length(flnm) == 62
    %         sub_name = flnm(end-5:end-4);
    %     elseif length(flnm) == 63
    %         sub_name = flnm(end-6:end-4);
    %     end
    
    % getting subject's name for old_hor
    %     switch role
    %         case 'predator'
    %             if length(flnm) == 72
    %                 sub_name = flnm(end-14:end-13);
    %             elseif length(flnm) == 71
    %                 sub_name = flnm(end-13);
    %             end
    %         case 'prey'
    %             if length(flnm) == 68
    %                 sub_name = flnm(end-10:end-9);
    %             elseif length(flnm) == 67
    %                 sub_name = flnm(end-9);
    %             end
    %     end
    
    %     % getting subject's name for ital
    %     if length(flnm) == 59
    %         sub_name = flnm(end-4);
    %     elseif length(flnm) == 60
    %         sub_name = flnm(end-5:end-4);
    %     end
    
    % getting subject's name for OT
    switch role
        case 'predator'
            k_pred   = k_pred+1;
            if length(flnm) == 72
                sub_name = flnm(end-18:end-13);
            end
            for k_model = 1:length(nmodel_array)
                for kk = 1:length(pred_PE(ii, :, k_model))
                    sub_tc                                  = sub_tc+1;
                    fitted_cell{sub_tc, 1}                  = sub_name;
                    fitted_cell{sub_tc, 2}                  = role;
                    fitted_cell{sub_tc, 3}                  = k_model;
                    fitted_cell{sub_tc, 4}                  = pred_EV(k_pred, kk, k_model);
                    fitted_cell{sub_tc, 5}                  = pred_PE(k_pred, kk, k_model);
                    fitted_cell{sub_tc, 6}                  = pred_risk(k_pred, kk, k_model);
                    fitted_cell{sub_tc, 7}                  = pred_risk_pe(k_pred, kk, k_model);
                end
            end
        case 'prey'
            k_prey   = k_prey+1;
            if length(flnm) == 68
                sub_name = flnm(end-14:end-9);
            end
            for k_model = 1:length(nmodel_array)
                for kk = 1:length(prey_PE(ii, :, k_model))
                    sub_tc                                  = sub_tc+1;
                    fitted_cell{sub_tc, 1}                  = sub_name;
                    fitted_cell{sub_tc, 2}                  = role;
                    fitted_cell{sub_tc, 3}                  = k_model;
                    fitted_cell{sub_tc, 4}                  = prey_EV(k_prey, kk, k_model);
                    fitted_cell{sub_tc, 5}                  = prey_PE(k_prey, kk, k_model);
                    fitted_cell{sub_tc, 6}                  = prey_risk(k_prey, kk, k_model);
                    fitted_cell{sub_tc, 7}                  = prey_risk_pe(k_prey, kk, k_model);
                end
            end
    end
end

% % 
% fid = fopen('/Users/michaelgiffin/Carsten PhD/OT_data/RL/data/OT_latent_params.txt', 'w');
% for ii = 1:length(fitted_cell)+1 % +1 for header
%     if ii == 1
%         fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\n', fitted_header{1,:});
%     else
%         fprintf(fid, '%s\t%s\t%d\t%f\t%f\t%f\t%f\n', fitted_cell{ii-1, :});
%     end
% end