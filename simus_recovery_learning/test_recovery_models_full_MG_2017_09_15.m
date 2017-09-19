clear
close all
clc

% what I need to do next is to play around with different number of
% sessions, and (first) different number of levels on the offers, and plot
% them all with the different recoveries and identifiabilites of each,
% including the R and P values, and save them somewhere, and then turn them
% into a powerpoint so that we can all have a look

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims  = 10; %40;                           % nsubs to simulates
n_trial = 10; %24;                           % ntrial per cond per session
n_sess  = 3;                            % nsession
offers  = 0:.5:10;
endow   = 10*ones(1,numel(offers));% parameters of the simulation

nmodel = 1;
% set up conditions and mutliple sessions
%------------------------------------------
Ra      = repmat(-[10,5,0],1,n_sess);            % predator true accepance thereshold (logit intercept)
Rb      = repmat([3,3,3],1,n_sess);              % predator true accpetance moise (logit slope)
n_cond  = size(Ra,2);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
Pa_rnd          = -10 + 10*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 1+2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
% Px_rnd          = random('gamma',1.2,5.0,n_sims,1);         %  Proposer  rating temperature
% Plr1_rnd        = .5*random('beta',1.1,1.1,n_sims,1);           %  Proposer  learning rate
% Plr2_rnd        = .5*random('beta',1.1,1.1,n_sims,1);           %  Proposer  learning rate

Px_rnd          = 1+4*rand(n_sims,1);         %  Proposer  rating temperature
Plr1_rnd        = .1+.8*rand(n_sims,1);           %  Proposer  learning rate
Plr2_rnd        = .1+.8*rand(n_sims,1);           %  Proposer  learning rate


% setup estimation
%---------------------
options         = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 1000000);
parameters      = NaN(n_sims,3,4,4);        
parametersLPP   = NaN(n_sims,3,4,4);
ll              = NaN(n_sims,4,4);        
LPP             = NaN(n_sims,4,4);

% Sim loop
for ktm = 1:4  % ktm = k true model
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
        
        n_rep           = 5;
        
        % for models 1 and 3 there is only 1 learning rate, for models 2
        % and 4 there are two learning rates
        switch ktm
            case 1 & 3
                [O,D] = learning_models_timeseries_1lr_MG_2017_09_18([bX,lr1],[Ra;Rb],n_trial,a0,b0,ktm);
                parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,2);
                lb = [0 0];
                LB = [0 0];
                ub = [15 1];
                UB = [Inf 1];
                num_params = 2;
            case 2 & 4
                [O,D] = learning_models_timeseries_2lr_MG_2017_09_18([bX,lr1,lr2],[Ra;Rb],n_trial,a0,b0,ktm);
                parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
                lb = [0 0 0];
                LB = [0 0 0];
                ub = [15 1 1];
                UB = [Inf 1 1];
                num_params = 3;
        end
       
        

        ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
        
%         lb = [0 0 0];          
%         LB = [0 0 0];
%         ub = [15 1 1];         
%         UB = [Inf 1 1];
        ddb = ub - lb;
        
        for kem = 1:4
            fprintf('executing loop 1, model %d of %d, sim %d of %d, and kem %d of %d\n', ktm, 4, k_sim, n_sims, kem, 4);
            
            for k_rep = 1:n_rep
%                 x0 = lb + rand(1,3).*ddb;
                switch kem
                    case 1 & 3
                        x0 = lb + rand(1,2).*ddb;
                        % %standard estimation
                        [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim_1lr_MG_2017_09_18(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                        % %lalace estimation
                        [parametersLPP_rep(k_rep,1:2),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning_MG_2017_09_18(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                    case 2 & 4
                        x0 = lb + rand(1,3).*ddb;
                        % %standard estimation
                        [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim_2lr_MG_2017_09_18(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                        % %lalace estimation
                        [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning_MG_2017_09_18(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                end
                x0 = lb + rand(1,num_params).*ddb;
                % %standard estimation
                [parameters_rep(k_rep,1:num_params),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                % %lalace estimation
                [parametersLPP_rep(k_rep,1:num_params),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
            end
            
            % here again we're having a problem with dimensionality,
            % because for some of the models we have 3 parameters and for
            % some we have 2, need to find an elegant way to fix this
            % (instead of the brute force tactic I've used above)
            [~,pos] = min(ll_rep);
            parameters(k_sim,:,ktm,kem)    =   parameters_rep(pos(1),:);
            ll(k_sim,ktm,kem)              =   ll_rep(pos(1),:);
            
            [~,posLPP] = min(LPP_rep);
            parametersLPP(k_sim,:,ktm,kem)      =   parametersLPP_rep(posLPP(1),:);
            LPP(k_sim,ktm,kem)                  =   LPP_rep(posLPP(1),:);
        end
    end
end

%%
for k_true = 1:4
    
    MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
    LL = squeeze(ll(:,k_true,:));
    
    nfpm=[2 3 2 3];
    
    for k_est=1:4
        bic(:,k_est)=-2*-LL(:,k_est)+nfpm(k_est)*log(3*n_trial*n_sess); % l2 is already positive
    end
    [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
    % [postBMC,outBMC]=VBA_groupBMC(-LL');
    BMC_output(k_true).post = postBMC;
    BMC_output(k_true).out = outBMC;
    
end


% save('test2')


%%
for k_true = 1:4
    
    model_titles = {'choice learning', 'choice learning 2 learning rates', 'reward learning', 'reward learning 2 learning rates'};
    legB = {'rating temperature','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    for k = 1:3
        
        subplot(2,3,k)
        plot(MP(:,k),squeeze(parameters(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        if k == 2
            title(model_titles{k_true})
        end
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(MP(:,k),parameters(:,k));
        
        subplot(2,3,3+k)
        plot(MP(:,k) ,squeeze(parametersLPP(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),parametersLPP(:,k));
        
    end
    
end