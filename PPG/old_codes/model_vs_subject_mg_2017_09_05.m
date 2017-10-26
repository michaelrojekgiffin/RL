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
nsub     = 5;   % number of plots to make, maximum is 166 (i.e. length(fl_dir) )
plot_ind = true;
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
    
    fprintf('estimating for subject %s, %s, %d out of %d\n', fl_dir(k_sub).name(6:11), predprey, k_sub, nsub);
    
    sub_o            = data(:,5);    % subject offer
    sub_r            = data(:,11);   % subject win/lose (logical)
    opponent_o       = data(:,6);    % column 6 is the choice of the opponent
    
    ntrial           = length(data);
    
    n_rep            = 10;
    parameters_rep   = NaN(n_rep,2);
    ll_rep           = NaN(n_rep,1);
    
    lb = [0 0];
    ub = [20 10];
    ddb = ub - lb;
    
    for k_rep = 1:n_rep
        x0 = lb + rand(1,2).*ddb;
        [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_2params_2017_09_04(x,priors.parameters(1:2),sub_o,sub_r,1,predprey),x0,[],[],[],[],lb,ub,[],options);
    end
    
    [~,pos] = min(ll_rep);
    
    parameters(k_sub,:)    =   parameters_rep(pos(1),:); % column 1 is beta (temperature), column 2 is alpha (learning rate)
    ll(k_sub)              =   ll_rep(pos(1),:);
    
    % Get slope and intercept of distribution against which the subject is
    % playind (i.e. get their true "acceptance function")
    opponent_o_freq = zeros(1, 11);
    opponent_probs = zeros(length(opponent_o_freq), 1);
    for k = 0:10
        opponent_o_freq(k+1)  = (sum(opponent_o==k)) / length(opponent_o);
        opponent_probs(k+1)   = sum(opponent_o_freq(1:k+1));
    end
    
    %     figure
    %     plot(opponent_o_freq);
    %
    % preallocate for simulation
    B0      = parameters(k_sub, 1);      % proposer temperature, estimated from subject data (specific to each subject)
    a0      = parameters(k_sub, 2);      % proposer learning rate, estimated from subject data (specific to each subject)
    v0      = priors.parameters(1);      % proposer initial prior on threshold (estimated from concatenated first trial investments)
    B       = priors.parameters(2);      % proposer noise (proposer's estimate of slope of opponent, (estimated from concatenated first trial investments))
    
    % pre-allocat
    %--------------------------
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
    V       = NaN(ntrial,1);             % selected offer
    kO      = NaN(ntrial,1);             % selected offer, in the 1:numel(offer) spavce
    O       = NaN(ntrial,1);             % selected offer, in the euro offer spavce
    Pd      = NaN(ntrial,1);             % trial  probability of accepting the selected offer
    D       = NaN(ntrial,1);             % final decision opf accepting the selected offer
    CPE     = NaN(ntrial,1);             % Choice Prediction error
    
    % initialize simulation
    %--------------------------
    V(1) = v0;
    for t = 1:ntrial
        
        % Proposer estimate the decision situation
        %-----------------------------------------------
        PA(t,:)             = logitp([V(t,:),B],offers);     % compute proba of accepting the offers given current model
        switch predprey
            case 'prey'
                EV(t,:)     = (endow - offers).* PA(t,:);                       % compute EV of the offers given current model
            case 'predator'
                % predator has adding value since they're guaranteed to
                % keep remainging endowment
                EV(t,:)     = (endow - offers) + ((endow - offers).* PA(t,:));  % compute EV of the offers given current model
        end
        
        % Proposer select an Offer
        %-----------------------------------------------
        % Soft max choice (using multinomial choice function)
        p       = exp(B0.*EV(t,:)) ./ sum(exp(B0.*EV(t,:)));      % multinomial choice function
        pd      = makedist('multinomial','probabilities',p);      % estimate the pdf from the pC
        ypdf    = pdf(pd,1:numel(offers));                        % generate pdf for the offers
        kO(t)   = random(pd);                                     % select index of pd based on probability created by makedist
        O(t)    = offers(kO(t));                                  % resample Offer in pdf (="soft-max")
        
        % D is whether the investment (proposal of our simulated subject)
        % is "accepted" (kill/survive) or "rejected" (miss/die)
        switch predprey
            case 'prey'
                D(t) = double(O(t) >= opponent_o(t));
            case 'predator'
                D(t) = double(O(t) > opponent_o(t));
        end
        
        
        % V is the "proposer"'s estimate of the "receiver's" intercept, and
        % B is the proposer's estimate of the receiver's slope. Pc is the
        % phi from powerpoint, it's the proposer's estimate of the
        % reciever's choice function, consisting of the intercept and
        % slope. Importantly, V (intercept) is updated on every trial with
        % rescorla-wagner update equation, while B (slope) remains constant
        Pc(t)       = logitp([V(t,:),B],O(t));
        
        % Updating Proposer estimation of the reciever's acceptance
        % function
        %------------------------------------------------------------
        CPE(t) = D(t) - Pc(t);
        V(t+1,:) = V(t,:) + a0.*CPE((t)); % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        
    end
    
    % record all parameters and all subject data for posterity
    switch predprey
        case 'predator'
            pred_ameters(k_pred, :)        = parameters(k_sub, :); % column 1 is temp, column 2 learning rate
            pred_dist(k_pred, :)           = sub_o;
            pred_sim_dist(k_pred, :)       = O;
            pred_opponent_probs(k_pred, :) = opponent_probs;
        case 'prey'
            prey_ameters(k_prey, :)        = parameters(k_sub, :); % column 1 is temp, column 2 learning rate
            prey_dist(k_prey, :)           = sub_o;
            prey_sim_dist(k_prey, :)       = O;
            prey_opponent_probs(k_prey, :) = opponent_probs;
    end
    
    % plot everything, per individual if plot_ind = true
    if plot_ind
        figure
        subplot(6, 3, 1:3:7)
        plot(1:length(data), sub_o, 'b-o', ...
            1:length(data), opponent_o, 'r--*',...
            1:length(data), O, '-.ks')
        legend(['Subject/' predprey], ['Opponent/' opponent], 'Simulated')
        title(['Subject (' predprey ') opponent (' opponent, ') and simiulation']);
        
        subplot(6, 3, 10:3:18)
        histogram(sub_o, 11)
        hold on
        histogram(O, 11)
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
            histogram(O(inv_counter:inv_row,:), 11)
            if ii == 1
                title('Real vs. simulated subject historgram, every 10 trials')
                legend('Subject', 'Simulation')
            end
            hold off
            inv_counter = inv_counter+10;
        end
        
        subplot(6, 3, 3:3:18)
        hold on
        plot(offers,opponent_probs,'-r')
        plot(offers,logitp(priors.parameters, offers),'-.k', 'Linewidth', 1.5)
        plot(offers,logitp([V(end), B(end)],offers),'-b')
        legend('true \phi','prior \phi','posterior \phi')
        xlabel('Offer (euro)')
        ylabel('p(Accept)')
        
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
