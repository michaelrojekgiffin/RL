% this script plots the "true" acceptance function of the opponent that the
% subject is playing against
%=========================================================================
%-------------------------------------------------------------------------
% EDIT - this part needs to be tweaked, depending on what I'm looking for
%-------------------------------------------------------------------------
nsub     = 10;   % number of plots to make, maximum is 166 (i.e. length(fl_dir) )

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%=========================================================================

close all

% define function
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));

options      = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);

cur_dir      = pwd;
data_dir     = fullfile(cur_dir,'data_matlab');
fl_dir       = dir(strcat(data_dir,filesep,'DATA_sub*'));

k_prey       = 0;
k_pred       = 0;

pred_ameters = [];
prey_ameters = [];
prey_freq    = [];
pred_freq    = [];
pred_dist    = [];
prey_dist    = [];
pred_sim_parameters = [];
prey_sim_parameters = [];


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
    sub_o            = data(:,5);    % subject offer
    sub_r            = data(:,11);   % subject win/lose (logical)
    opponent_o       = data(:,6);    % column 6 is the choice of the opponent
    
    ntrial           = length(data);
    
%     n_rep            = 10;
%     parameters_rep   = NaN(n_rep,2);
%     ll_rep           = NaN(n_rep,1);
%     
%     lb = [0 0];
%     ub = [20 10];
%     ddb = ub - lb;
%     
%     for k_rep = 1:n_rep
%         x0 = lb + rand(1,2).*ddb;
%         [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_2params_2017_09_04(x,priors.parameters(1:2),sub_o,sub_r,1,predprey),x0,[],[],[],[],lb,ub,[],options);
%     end
%     
%     [~,pos] = min(ll_rep);
    
%     parameters(k_sub,:)    =   parameters_rep(pos(1),:); % column 1 is beta (temperature), column 2 is alpha (learning rate)
%     ll(k_sub)              =   ll_rep(pos(1),:);
%     
    % Get slope and intercept of distribution against which the subject is
    % playind (i.e. get their true "acceptance function")
    opponent_parameters    =    empirical_priors_mg_simVsub_05_09_2017(opponent_o, opponent);
    figure
    hold on
    plot(offers,logitp(opponent_parameters, offers),'-r')
    title(sprintf('%s of subject %d, slope: %.2f, intercept: %.2f', opponent, k_sub, opponent_parameters(1), opponent_parameters(2)));
    
end