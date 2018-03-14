% this script fits the model to subject data, and then runs a simulation
% for each subject using the parameters estimated from their data in order
% to compare
clear all
close all
clc
%=========================================================================
%-------------------------------------------------------------------------
% EDIT - this part needs to be tweaked, depending on what I'm looking for
%-------------------------------------------------------------------------
nsub          = 1;   % number of plots to make, maximum is 166 (i.e. length(fl_dir) )
plot_ind      = true;
nsim          = 15;
plot_all_data = false;
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%=========================================================================

% define function
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));

options      = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000);

cur_dir      = pwd;
data_dir     = fullfile(cur_dir,'data_matlab');
fl_dir       = dir(strcat(data_dir,filesep,'DATA_sub*'));

k_prey       = 0;
k_pred       = 0;

pred_ameters        = []; % parameters of predator model for each subject playing as predator
prey_ameters        = []; % parameters of prey model for each subject playing as prey
prey_freq           = []; % prey frequencies
pred_freq           = []; % predator frequencies
pred_dist           = []; % predator distribution
prey_dist           = []; % prey distribution
pred_sim_dist       = []; % distribution created from the simulations of predator (for each subject)
prey_sim_dist       = []; % distribution created from the simulations of prey (for each subject)
pred_opponent_probs = []; % distribution each predator played against
prey_opponent_dist  = []; % distribution each prey played against

EV_sub              = NaN(nsub, 60);             % expected value of each decision given parameters recovered from data (each subject has own row, each trial has own column)
PA_sub              = NaN(nsub, length(offers)); % probaility that each investment will succeed given params of model
V_sub               = NaN(nsub, 1);              % posterior intercept of subject decision function given params of model
EV_sub_posterior    = NaN(nsub, length(offers)); % expected value of every possible offer with posterior intercept given recovered model params
for k_sub    = 1:nsub
    
    flnm     = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
            k_pred   = k_pred+1;
            predprey = 'predator';
            opponent = 'prey';
            priors   = load('data_matlab/predator_priors.mat');
        case 'prey'
            k_prey   = k_prey+1;
            predprey = 'prey';
            opponent = 'predator';
            priors   = load('data_matlab/prey_priors.mat');
    end
%     
    fprintf('estimating for subject %s, %s, %d out of %d\n', fl_dir(k_sub).name(6:11), predprey, k_sub, nsub);
    
    sub_o            = data(:,5);    % subject offer
    sub_r            = data(:,11);   % subject win/lose (logical)
    opponent_o       = data(:,6);    % column 6 is the choice of the opponent
    
    ntrial           = length(data);
    
    n_rep            = 50;
    parameters_rep   = NaN(n_rep,2);
    ll_rep           = NaN(n_rep,1);
    
    lb = [0 0];
    ub = [20 11];
    ddb = ub - lb;
    
    for k_rep = 1:n_rep
        x0 = lb + rand(1,2).*ddb;
        [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_2params_2017_09_11(x,priors.parameters(1:2),sub_o,sub_r,1,predprey),x0,[],[],[],[],lb,ub,[],options);
    end
    
    [~,pos] = min(ll_rep);
    
    parameters(k_sub,:)    =   parameters_rep(pos(1),:); % column 1 is beta (temperature), column 2 is alpha (learning rate)
    ll(k_sub)              =   ll_rep(pos(1),:);
    
    % record subject's expected values (sub_EV), estimates of probability
    % of success (sub_PA), and intercept (sub_V), given the parameters
    % estimated above
    [~, all_EV, all_PA, all_V] = ModelEstimation_2params_2017_09_11(parameters(k_sub,:), priors.parameters(1:2),sub_o,sub_r,1,predprey);
    V_sub(k_sub)               = all_V(end);     % posterior intercept
    PA_sub(k_sub, :)           = all_PA(end, :); % posterior prob of success
    EV_sub_posterior(k_sub, :) = all_EV(end, :); % posterior EV of all offers
    
    for ee = 1:length(all_EV)
        EV_sub(k_sub, ee) = all_EV(ee, sub_o(ee)+1); % EV of each offer selected by subject
    end
    
    
    % Get slope and intercept of distribution against which the subject is
    % playind (i.e. get their true "acceptance function")
    opponent_o_freq = zeros(1, 11);
    succes_probs = zeros(length(opponent_o_freq), 1); % probability of success
    for k = 0:10
        opponent_o_freq(k+1)       = (sum(opponent_o==k)) / length(opponent_o);
        switch predprey
            case 'prey'
                succes_probs(k+1) = sum(opponent_o_freq(1:k+1));
            case 'predator'                                         % predator only wins if prey invests less
                succes_probs(k+1) = sum(opponent_o_freq(1:k));     % sum all probabilities below chosen offer
        end
    end
    empirical_EV = zeros(length(offers), 1);
    switch predprey
        case 'prey'
            empirical_EV = (endow - offers).* succes_probs';
        case 'predator'
            empirical_EV = (endow - offers) + ((endow - offers).* succes_probs');
    end
    
    % preallocate for simulation
    B0      = parameters(k_sub, 1);      % proposer temperature, estimated from subject data (specific to each subject)
    a0      = parameters(k_sub, 2);      % proposer learning rate, estimated from subject data (specific to each subject)
    v0      = priors.parameters(1);      % proposer initial prior on threshold (estimated from concatenated first trial investments)
    B       = priors.parameters(2);      % proposer noise (proposer's estimate of slope of opponent, (estimated from concatenated first trial investments))
    
    % pre-allocat,
    %--------------------------
    PA      = NaN(ntrial,numel(offers), nsim); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers), nsim); % estimated expected value of all offers
    V       = NaN(ntrial, nsim);               % estimated opponent intercept
    kO      = NaN(ntrial, nsim);               % selected offer, in the 1:numel(offer) space
    O       = NaN(ntrial, nsim);               % selected offer, in the euro offer spavce
    Pd      = NaN(ntrial, nsim);               % trial  probability of accepting the selected offer
    D       = NaN(ntrial, nsim);               % final decision of accepting the selected offer
    CPE     = NaN(ntrial, nsim);               % Choice Prediction error
    
    % initialize simulation
    %--------------------------
    for i_sim = 1:nsim
        V(1, i_sim) = v0;
        for t = 1:ntrial
            
            % Proposer estimate the decision situation
            %-----------------------------------------------
            PA(t,:, i_sim)             = logitp([V(t,i_sim),B],offers);     % compute proba of accepting the offers given current model
            switch predprey
                case 'prey'
                    EV(t,:, i_sim)     = (endow - offers).* PA(t,:, i_sim); % compute EV of the offers given current model
                case 'predator'
                    % predator has added value since they're guaranteed to
                    % keep remainging endowment
                    EV(t,:, i_sim)     = (endow - offers) + ((endow - offers).* PA(t,:, i_sim));  % compute EV of the offers given current model
            end
            
            % Proposer select an Offer
            %-----------------------------------------------
            % Soft max choice (using multinomial choice function)
            p              = exp(B0.*EV(t,:, i_sim)) ./ sum(exp(B0.*EV(t,:, i_sim)));  % multinomial choice function
            pd             = makedist('multinomial','probabilities',p);                % estimate the pdf from the pC
            ypdf           = pdf(pd,1:numel(offers));                                  % generate pdf for the offers
            kO(t, i_sim)   = random(pd);                                               % select index of pd based on probability created by makedist
            O(t, i_sim)    = offers(kO(t, i_sim));                                     % resample Offer in pdf (="soft-max")
            
            % D is whether the investment (proposal of our simulated subject)
            % is "accepted" (kill/survive) or "rejected" (miss/die)
            switch predprey
                case 'prey'
                    D(t, i_sim) = double(O(t, i_sim) >= opponent_o(t));
                case 'predator'
                    D(t, i_sim) = double(O(t, i_sim) > opponent_o(t));
            end
            
            % V is the "proposer"'s estimate of the "receiver's" intercept, and
            % B is the proposer's estimate of the receiver's slope. Pc is the
            % phi from powerpoint, it's the proposer's estimate of the
            % reciever's choice function, consisting of the intercept and
            % slope. Importantly, V (intercept) is updated on every trial with
            % rescorla-wagner update equation, while B (slope) remains constant
            Pc(t, i_sim)       = logitp([V(t, i_sim),B],O(t, i_sim));
            
            % Updating Proposer estimation of the reciever's acceptance
            % function
            %------------------------------------------------------------
            CPE(t, i_sim) = D(t, i_sim) - Pc(t, i_sim);
            V(t+1,i_sim) = V(t,i_sim) + a0.*CPE((t), i_sim); % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        end
    end
    % average across all simulations
    V_sim              = mean(V(end, :)); % average of V (opponent's intercept) estimate across simulations
    % pre-allocate
    EV_sim_prior       = NaN(1,length(offers));
    EV_sim_posterior   = NaN(1,length(offers));
    PA_sim_prior       = NaN(1,length(offers));
    PA_sim_posterior   = NaN(1,length(offers));
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
        O_mean(m_trial) = mean(O(m_trial, :));
    end
    % record all parameters and all subject data for posterity
    switch predprey
        case 'predator'
            pred_ameters(k_pred, :)        = parameters(k_sub, :); % column 1 is temperature, column 2 learning rate
            pred_dist(k_pred, :)           = sub_o;
            pred_sim_dist(k_pred, :)       = O_mean;
            pred_opponent_probs(k_pred, :) = succes_probs;
        case 'prey'
            prey_ameters(k_prey, :)        = parameters(k_sub, :); % column 1 is temperature, column 2 learning rate
            prey_dist(k_prey, :)           = sub_o;
            prey_sim_dist(k_prey, :)       = O_mean;
            prey_opponent_probs(k_prey, :) = succes_probs;
    end
    
    % plot everything, per individual if plot_ind = true
    if plot_ind
        if ~plot_all_data
            figure
            subplot(1, 2, 1)
            hold on
            plot(offers,succes_probs,'-r')
            plot(offers,logitp(priors.parameters(1:2), offers),'-.k', 'Linewidth', 1.5)
            plot(offers,logitp([V_sub(k_sub), B(end)],offers),'--b')
            plot(offers,logitp([V_sim, B],offers),'--g')
            legend('true \phi','prior \phi','posterior \phi', 'simulated posterior \phi')
            xlabel('Offer (euro)')
            ylabel('p(Accept)')
            title(['Subject number ' k_sub ' (' predprey ')']);
            
            % this plot is wrong, it is depicting simulation, should be
            % depciting subject - need to simply add a line for the subject
            subplot(1, 2, 2)
            hold on
            plot(offers, empirical_EV, '-r')
            plot(offers, EV_sim_prior, '-.k', 'Linewidth', 1.5)
            plot(offers, EV_sub_posterior(k_sub, :), '--b')
            plot(offers, EV_sim_posterior, '--g')
            legend('true', 'prior', 'posterior', 'simulated posterior')
            xlabel('Offer (euro)')
            ylabel('Expected Value')
            
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
            plot(offers,succes_probs,'-r')
            plot(offers,logitp(priors.parameters(1:2), offers),'-.k', 'Linewidth', 1.5)
            plot(offers,logitp([V_sub(k_sub), B(end)],offers),'--b')
            plot(offers,logitp([V_sim, B],offers),'--g')
            legend('true \phi','prior \phi','posterior \phi', 'simulated posterior \phi')
            xlabel('Offer (euro)')
            ylabel('p(Accept)')
            title(['Subject (' predprey ')']);
            
            subplot(6, 3, 12:3:18)
            hold on
            plot(offers, empirical_EV, '-r')
            plot(offers, EV_sim_prior, '-.k', 'Linewidth', 1.5)
            plot(offers, EV_sub_posterior(k_sub, :), '--b')
            plot(offers, EV_sim_posterior, '--g')
            legend('true', 'prior', 'posterior', 'simulated posterior')
            xlabel('Offer (euro)')
            ylabel('Expected Value')
            
        end
    end
    
end

%% write to file
% write to file for analysis in R
% [H P] = ttest(pred_ameters(2, :), prey_ameters(2, :))
% [H P] = ttest(pred_ameters(1, :), prey_ameters(1, :))

% fitted_header = {'sub_name', 'alpha', 'beta'};
% fitted_length = nsub;
% fitted_cell = cell(fitted_length, length(fitted_header));
%
% % % % sub_counter = 1;
% % % % sub_tc = 0; % subject trial counter
% k_pred = 0;
% k_prey = 0;
% for ii = 1:fitted_length
%     sub_tc   = sub_tc+1;
%     flnm     = fullfile(data_dir,fl_dir(ii).name);
%     load(flnm)
%     if length(flnm) == 57
%         sub_name = flnm(end-4);
%     elseif length(flnm) == 58
%         sub_name = flnm(end-5:end-4);
%     elseif length(flnm) == 59
%         sub_name = flnm(end-6:end-4);
%     end
%
%     switch role
%         case 'predator'
%             k_pred   = k_pred+1;
%             fitted_cell{ii, 1} = sub_name; % column 2 of pred_ameters is alpha
%             fitted_cell{ii, 2} = pred_ameters(k_pred, 2); % column 2 of pred_ameters is alpha
%             fitted_cell{ii, 3} = pred_ameters(k_pred, 1); % column 1 of pred_ameters is beta
%         case 'prey'
%             k_prey   = k_prey+1;
%             fitted_cell{ii, 1} = sub_name; % column 2 of pred_ameters is alpha
%             fitted_cell{ii, 2} = prey_ameters(k_prey, 2); % column 2 of pred_ameters is alpha
%             fitted_cell{ii, 3} = prey_ameters(k_prey, 1); % column 1 of pred_ameters is beta
%     end
% end
%
%
%
% fid = fopen('/Users/michaelgiffin/Carsten PhD/hormones/data/modeling/hor_fitted.txt', 'w');
% for ii = 1:length(fitted_cell)+1 % +1 for header
%     if ii == 1
%         fprintf(fid, '%s\t%s\t%s\n', fitted_header{1,:});
%     else
%         fprintf(fid, '%s\t%f\t%f\n', fitted_cell{ii-1, :});
%     end
% end
%
