% this makes some plots to explore whats going on over time in the hor
% dataset. 

% after this what I'd like to do is to first see what sort of distribution
% my parameters exhibit, and then fit the models again, and even better
% would be to fit one model to both predator and prey, varying the
% different parameters between the two of them, and see which model
% explains the data the best - then I can directly compare the parameters
% between conditions because it will be the same model

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
ste = @(x) (std(x))/(sqrt(length(x)));
offers = 0:1:10;
endow  = 10*ones(1,numel(offers));% parameters of the simulation

ntr = 60;
% pre-allocate
pred_offer_mat      = NaN(nsub, ntr);
pred_accept_mat     = NaN(nsub, ntr);
pred_reward_mat     = NaN(nsub, ntr);
pred_RT_mat         = NaN(nsub, ntr);

prey_offer_mat      = NaN(nsub, ntr);
prey_accept_mat     = NaN(nsub, ntr);
prey_reward_mat     = NaN(nsub, ntr);
prey_RT_mat         = NaN(nsub, ntr);
for k_sub = 1:nsub
    
    flnm = fullfile(data_dir,fl_list(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
            pred_offer_mat(k_sub, :)   =  data(:, 5)';  % column 5 is playerselect
            pred_accept_mat(k_sub, :)  =  data(:, 11)'; % column 11 is payoff
            pred_reward_mat(k_sub, :)  =  data(:, 9)';  % column 9 is payoff
            pred_opponent_mat(k_sub, :) = data(:, 6)'; % the investment of the opponent
            
            pred_RT_mat(k_sub, :)       = data(:, 7)'; % the RT
            
            O                          = squeeze(pred_offer_mat(k_sub,:));
            D                          = squeeze(pred_accept_mat(k_sub,:));
            R_o                        = data(:, 6)';   % column 6 is offer of opponent
            
            role_dim                   = 1;             % tells which dimension to store loglik and parameters
        case 'prey'
            prey_offer_mat(k_sub, :)   =  data(:, 5)';  % column 5 is playerselect
            prey_accept_mat(k_sub, :)  =  data(:, 11)'; % column 11 is payoff
            prey_reward_mat(k_sub, :)  =  data(:, 9)';  % column 9 is payoff
            prey_opponent_mat(k_sub, :) = data(:, 6)'; % the investment of the opponent
            prey_RT_mat(k_sub, :)       = data(:, 7)'; % the RT
            
            O                          = squeeze(prey_offer_mat(k_sub,:));
            D                          = squeeze(prey_accept_mat(k_sub,:));
            R_o                        = data(:, 6)';   % column 6 is offer of opponent
            
            role_dim                   = 2;             % tells which dimension to store loglik and parameters
    end
    
   
    
    
end
%% plotting investments
% ----------------------------------------------------------
% plotting investments
% ----------------------------------------------------------
close all
pred_offer_mat(any(isnan(pred_offer_mat), 2), :) = [];
pred_opponent_mat(any(isnan(pred_opponent_mat), 2), :) = [];


prey_offer_mat(any(isnan(prey_offer_mat), 2), :) = [];
prey_opponent_mat(any(isnan(prey_opponent_mat), 2), :) = [];

count = 1;
figure 
subplot(1, 3, 1:2)
hold on
for ii = 1:size(pred_offer_mat, 2)
%     plot(ii, mean(pred_offer_mat(:, ii)),  'r--o', 'MarkerFaceColor','red')
%     plot(ii, mean(prey_offer_mat(:, ii)),  'b--d', 'MarkerFaceColor','blue')
    
    plot(ii, mean(pred_offer_mat(:, ii)),  'r--o')
    plot(ii, mean(prey_offer_mat(:, ii)),  'b--d')
    if mod(ii, 10) == 0
        errorbar(ii-6, mean(mean(pred_offer_mat(:, count:count+9))), mean(ste(pred_offer_mat(:, count:count+9))),...
            'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
        errorbar(ii-4, mean(mean(prey_offer_mat(:, count:count+9))), mean(ste(prey_offer_mat(:, count:count+9))),...
        'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
        count = count + 10;
    end
end
hold off


subplot(1, 3, 3)
hold on
errorbar(1, mean(mean(pred_offer_mat(:, :))), mean(ste(pred_offer_mat(:, :))),...
    'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
errorbar(2, mean(mean(prey_offer_mat(:, :))), mean(ste(prey_offer_mat(:, :))),...
    'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
xlim([.8, 2.2])
hold off

%% rewards
% ----------------------------------------------------------
% plotting rewards
% ----------------------------------------------------------

pred_reward_mat(any(isnan(pred_reward_mat), 2), :) = [];
pred_opponent_mat(any(isnan(pred_opponent_mat), 2), :) = [];


prey_reward_mat(any(isnan(prey_reward_mat), 2), :) = [];
prey_opponent_mat(any(isnan(prey_opponent_mat), 2), :) = [];

count = 1;
figure 
subplot(1, 3, 1:2)
hold on
for ii = 1:size(pred_reward_mat, 2)
%     plot(ii, mean(pred_reward_mat(:, ii)),  'r--o', 'MarkerFaceColor','red')
%     plot(ii, mean(prey_reward_mat(:, ii)),  'b--d', 'MarkerFaceColor','blue')
    
    plot(ii, mean(pred_reward_mat(:, ii)),  'r--o')
    plot(ii, mean(prey_reward_mat(:, ii)),  'b--d')
    if mod(ii, 10) == 0
        errorbar(ii-6, mean(mean(pred_reward_mat(:, count:count+9))), mean(ste(pred_reward_mat(:, count:count+9))),...
            'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
        errorbar(ii-4, mean(mean(prey_reward_mat(:, count:count+9))), mean(ste(prey_reward_mat(:, count:count+9))),...
        'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
        count = count + 10;
    end
end
hold off


subplot(1, 3, 3)
hold on
errorbar(1, mean(mean(pred_reward_mat(:, :))), mean(ste(pred_reward_mat(:, :))),...
    'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
errorbar(2, mean(mean(prey_reward_mat(:, :))), mean(ste(prey_reward_mat(:, :))),...
    'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
xlim([.8, 2.2])
hold off


%% RT
% ----------------------------------------------------------
% plotting rewards
% ----------------------------------------------------------

pred_RT_mat(any(isnan(pred_RT_mat), 2), :) = [];
pred_opponent_mat(any(isnan(pred_opponent_mat), 2), :) = [];


prey_RT_mat(any(isnan(prey_RT_mat), 2), :) = [];
prey_opponent_mat(any(isnan(prey_opponent_mat), 2), :) = [];

count = 1;
figure 
subplot(1, 3, 1:2)
hold on
for ii = 1:size(pred_RT_mat, 2)
%     plot(ii, mean(pred_RT_mat(:, ii)),  'r--o', 'MarkerFaceColor','red')
%     plot(ii, mean(prey_RT_mat(:, ii)),  'b--d', 'MarkerFaceColor','blue')
    
    plot(ii, mean(pred_RT_mat(:, ii)),  'r--o')
    plot(ii, mean(prey_RT_mat(:, ii)),  'b--d')
    if mod(ii, 10) == 0
        errorbar(ii-6, mean(mean(pred_RT_mat(:, count:count+9))), mean(ste(pred_RT_mat(:, count:count+9))),...
            'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
        errorbar(ii-4, mean(mean(prey_RT_mat(:, count:count+9))), mean(ste(prey_RT_mat(:, count:count+9))),...
        'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
        count = count + 10;
    end
end
ylim([0 8000])
hold off


subplot(1, 3, 3)
hold on
errorbar(1, mean(mean(pred_RT_mat(:, :))), mean(ste(pred_RT_mat(:, :))),...
    'r--o', 'MarkerSize',10, 'MarkerEdgeColor','red','MarkerFaceColor','red')
errorbar(2, mean(mean(prey_RT_mat(:, :))), mean(ste(prey_RT_mat(:, :))),...
    'b--d', 'MarkerSize',10, 'MarkerEdgeColor','blue','MarkerFaceColor','blue')
xlim([.8, 2.2])
ylim([1800 3000])
hold off
