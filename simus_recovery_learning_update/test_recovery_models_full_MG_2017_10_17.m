% clear
% close all
% clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims  = 20;                           % nsubs to simulates
% n_trial = 10;                           % ntrial per cond per session
n_sess          = 2;                            % nsession
offers          = 0:1:10;
endow           = 10*ones(1,numel(offers));% parameters of the simulation
lr_upper_bound  = 1;

trial_space = 30:10:60;

% Generate params
%-------------------
Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 2.9+.2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)

Px_rnd          = .5+2.5*rand(n_sims,1);         %  Proposer  rating temperature
% Px_rnd          = 3+3*rand(n_sims,1);         %  Proposer  rating temperature
Plr1_rnd        = rand(n_sims,lr_upper_bound);           %  Proposer  learning rate
Plr2_rnd        = rand(n_sims,1);           %  Proposer  learning rate



for kk_trial= 1:length(trial_space)
    n_trial = trial_space(kk_trial);
    
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
    
%     % Generate params
%     %-------------------
%     Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
%     Pb_rnd          = 2.9+.2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
%     
%     Px_rnd          = .5+2.5*rand(n_sims,1);         %  Proposer  rating temperature
%     % Px_rnd          = 3+3*rand(n_sims,1);         %  Proposer  rating temperature
%     Plr1_rnd        = rand(n_sims,1);           %  Proposer  learning rate
%     Plr2_rnd        = rand(n_sims,1);           %  Proposer  learning rate
%     
%     
    % setup estimation
    %---------------------
    options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);
    parameters  = NaN(n_sims,3,nmods,nmods);        parametersLPP  = NaN(n_sims,3,nmods,nmods);
    ll          = NaN(n_sims,nmods,nmods);        LPP            = NaN(n_sims,nmods,nmods);
    
   
    PE          = NaN(n_sims, n_trial, nmods); % prediction errors of models
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
            
            % simulates how subjects learn multiple conditions (multiple
            % resistance points) of responders given the number of trials we
            % stipulate
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
                    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim_MG(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                    % %lalace estimation
                    [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2_MG(x,O,D,a0,b0,kem, lr_upper_bound),x0,[],[],[],[],LB,UB,[],options);
                end
                [~,pos] = min(ll_rep);
                parameters(k_sim,:,ktm,kem)    =   parameters_rep(pos(1),:);
                ll(k_sim,ktm,kem)              =   ll_rep(pos(1),:);
                
                [~,posLPP] = min(LPP_rep);
                parametersLPP(k_sim,:,ktm,kem)      =   parametersLPP_rep(posLPP(1),:);
                LPP(k_sim,ktm,kem)                  =   LPP_rep(posLPP(1),:);
                
            end
            
            [~, PE(k_sim, :, ktm)]                 = learning_models_estim_MG(parametersLPP(k_sim,:,ktm,ktm),O,D,a0,b0,nmodel);
        end
    end
    
    % on each iteration, this calculates the BIC of the true model verses the
    % other models
    for k_true = modelspace
    
%         VBA_fig = figure;
        
        MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
    %     LL = squeeze(ll(:,k_true,:));
        LL = squeeze(LPP(:,k_true,:));
    
        for k_est= modelspace
            bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
        end
        [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
        BMC_output(k_true).post = postBMC;
        BMC_output(k_true).out = outBMC;
        
%         print(VBA_groupBMC(-bic'./2), ['..' filesep 'reports' filesep 'figures' filesep 'pres_2017_10_18' filesep, ...
%         'VBA_model_', num2str(k_true), num2str(n_sims), 'sims', num2str(n_trial), 'trials', ...
%         num2str(n_sess), 'sess', num2str(n_cond), 'cond' ], '-dpng');
    
    end
    
    
    BIC      = [];
    BIC_mean = [];
    for k_true = 1:4
        for k_est=1:4
            [~, BIC(:,k_true, k_est)] = aicbic(LPP(:, k_true, k_est), nfpm(k_est), n_trial);
            
            BIC_mean(k_true, k_est)          = mean(BIC(:, k_true, k_est));
            ll_mean(k_true, k_est)           = mean(-LPP(:, k_true, k_est));
        end
    end
    BIC_fig = figure;
    imagesc(BIC_mean)
    colsize = size(BIC_mean, 2);
    colcounter = 0;
    rowcounter = 1;
    for ii = 1:numel(BIC_mean)
        colcounter = colcounter+1;
        text(colcounter, rowcounter, num2str(BIC_mean(rowcounter, colcounter)), 'Color', 'r')
        if colcounter == colsize
            colcounter = 0;
            rowcounter = rowcounter + 1;
        end
    end
    colormap default
    title(sprintf('BIC'))
    colorbar
    print(BIC_fig, ['..' filesep 'reports' filesep 'figures' filesep 'pres_2017_10_18' filesep, ...
        'BIC_' num2str(n_sims), 'sims', num2str(n_trial), 'trials', ...
        num2str(n_sess), 'sess', num2str(n_cond), 'cond' ], '-dpng');
    
    
    
    for k_true = modelspace
        
        num_free_params = nfpm(k_true);
        legB            = {'rating temperature','learning rate 1','learning rate 2'};
        
        model_rec_fig = figure;
        set(gcf,'Color',[1,1,1])
        for k = 1:num_free_params
            
            %         subplot(2,num_free_params,k)
            %         plot(MP(:,k),squeeze(parameters(:,k,k_true,k_true)),'o',...
            %             'MarkerEdgeColor',[0,0,0],...
            %             'MarkerFaceColor',[1,1,1])
            %         xlabel(strcat(['true ' legB{k}]));
            %         ylabel(strcat(['estimated ' legB{k}]));
            %         lsline;
            %         [corrR(k),corrP(k)] = corr(MP(:,k),squeeze(parameters(:,k,k_true,k_true)));
            %
            %         txt1 = sprintf('r = %f\n p = %f', corrR(k),corrP(k));
            %         if k == 1
            %             title(sprintf('model %d, %d trials\n %s', k_true,n_trial, txt1));
            %         else
            %             title(txt1)
            %         end
            %
            
            subplot(1,num_free_params,k)
            plot(MP(:,k) ,squeeze(parametersLPP(:,k,k_true,k_true)),'o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,1])
            xlabel(strcat(['true ' legB{k}]));
            ylabel(strcat(['estimated ' legB{k} ' LPP']));
            lsline;
            
            [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),squeeze(parametersLPP(:,k,k_true,k_true)));
            txt1 = sprintf('r = %f\n p = %f', corrR_LPP(k),corrP_LPP(k));
            if k == 1
                title(sprintf('model %d, %d trials\n %s', k_true,n_trial, txt1));
            else
                title(txt1)
            end
            
        end
        
        print(model_rec_fig, ['..' filesep 'reports' filesep 'figures' filesep 'pres_2017_10_18' filesep 'model', ...
            num2str(k_true) 'recovery_' num2str(n_sims), 'sims', num2str(n_trial), 'trials', ...
            num2str(n_sess), 'sess', num2str(n_cond), 'cond' ], '-dpng');
        
    end
    
end
% 
% PE1 = reshape(PE(:, :, 1), [],1);
% PE2 = reshape(PE(:, :, 2), [],1);
% PE3 = reshape(PE(:, :, 3), [],1);
% PE4 = reshape(PE(:, :, 4), [],1);