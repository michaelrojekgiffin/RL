clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims      = 40;                           % nsubs to simulates
n_trial     = 18;                           % ntrial per cond per session
n_sess      = 4;                            % nsession
offers      = 0:1:10;
suboffers   = 0:.1:10;
endow       = 10*ones(1,numel(offers));% parameters of the simulation
subendow    = 10*ones(1,numel(suboffers));% parameters of the simulation
nmodel_array= 1:4; % all the models to loop through
% nmodel      = 2;
lr_upper_bound = 3;

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
% logitp = @(b,x) b(3)/2+(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));
%logitp = @(b,x)(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));
% Ubch = .10;

% Generate params
%-------------------
Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 2.5+1*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)

Px_rnd          = .5+2.5*rand(n_sims,1);         %  Proposer  rating temperature
% Px_rnd          = 3+3*rand(n_sims,1);         %  Proposer  rating temperature

Plr2_rnd        = rand(n_sims,1);           %  Proposer  learning rate
Plr1_rnd        = 3*rand(n_sims,1);           %  Proposer  learning rate

% PreAllocate
%---------------
A_mat = zeros(n_sims,n_trial,nc);
PE_mat = zeros(n_sims,n_trial,nc);
PAsub_mat = NaN(n_sims,numel(suboffers));

% setup estimation
%---------------------
options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000, 'display', 'off');
parameters  = NaN(n_sims,3);        parametersLPP  = NaN(n_sims,3);
ll          = NaN(n_sims,1);        LPP            = NaN(n_sims,1);

% for confusion matrix
% dim 1: simulated subject (so each entry is one simulated subject)
% dim 2: paraemters - so each entry is a different parameter of the model
% dim 3: the actual model we used to generate the simulated data
% dim 4: the parameters estimated from
con_ll              = NaN(n_sims,4,4,4);
con_LPP             = NaN(n_sims,4,4,4);

for model_iter = 1:length(nmodel_array)
    nmodel = nmodel_array(model_iter);
    % Sim loop
    %----------
    for k_sim = 1:n_sims
        
        fprintf('model %d, sim %d of %d\n', nmodel, k_sim, n_sims);
        
        % pre-allocate
        O_mat = NaN(n_trial,n_cond);
        D_mat = NaN(n_trial,n_cond);
        
        % get params
        a0  = Pa_rnd(k_sim);
        b0  = Pb_rnd(k_sim);
        bX  = Px_rnd(k_sim);
        lr1 = Plr1_rnd(k_sim);
        lr2 = Plr2_rnd(k_sim);
        
        [O,D,PE,at] = learning_models_timeseries_MG_2017_10_12([bX,lr1,lr2],[Ra;Rb],n_trial,a0,b0,nmodel);
        
        for k_sess = 1:n_sess
            k_in = (k_sess-1)*nc+1;
            k_out = k_sess*nc;
            A_mat(k_sim,:,:) = squeeze(A_mat(k_sim,:,:)) + at(1:n_trial,k_in:k_out)./n_sess;
            PE_mat(k_sim,:,:) = squeeze(PE_mat(k_sim,:,:)) + PE(1:n_trial,k_in:k_out)./n_sess;
        end
        
        PAsub_mat(k_sim,:) = logitp([a0,b0],suboffers);
        
        n_rep           = 5;
        parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
        ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
        
        lb = [0 0 0];          LB = [0 0 0];
        ub = [5 1 1];         UB = [Inf 3 1];
        ddb = ub - lb;
        
        for k_rep = 1:n_rep
            x0 = lb + rand(1,3).*ddb;
            %standard estimation
            [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,a0,b0,nmodel),x0,[],[],[],[],LB,UB,[],options);
            %lalace estimation
            [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2_MG(x,O,D,a0,b0,nmodel,lr_upper_bound),x0,[],[],[],[],LB,UB,[],options);
        end
        
        [~,pos] = min(ll_rep);
        parameters(k_sim,:)    =   parameters_rep(pos(1),:);
        ll(k_sim)              =   ll_rep(pos(1),:);
        
        [~,posLPP] = min(LPP_rep);
        parametersLPP(k_sim,:)      =   parametersLPP_rep(posLPP(1),:);
        LPP(k_sim)                  =   LPP_rep(posLPP(1),:);
        
        
        %%%%% for confusion matrices
        for con_model_iter = 1:length(nmodel_array)
            con_nmodel = nmodel_array(con_model_iter);
            
            if con_nmodel == 1 || con_nmodel == 3
                numfreeparams = 2;
                ub = [5 lr_upper_bound];         UB = [Inf lr_upper_bound];
            else
                numfreeparams = 3;
                ub = [5 lr_upper_bound 1];         UB = [Inf lr_upper_bound 1];
            end
            
            parameters_rep  = NaN(n_rep,numfreeparams);     parametersLPP_rep  = NaN(n_rep,numfreeparams);
            ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
            
            lb = zeros(1, numfreeparams);          LB = zeros(1, numfreeparams);
            
            ddb = ub - lb;
            
            for k_rep = 1:n_rep
                x0 = lb + rand(1,numfreeparams).*ddb;
                %standard estimation
                %             [parameters_rep(k_rep,1:numfreeparams),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,a0,b0,nmodel),x0,[],[],[],[],LB,UB,[],options);
                %lalace estimation
                [parametersLPP_rep(k_rep,1:numfreeparams),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2_MG(x,O,D,a0,b0,con_nmodel,lr_upper_bound),x0,[],[],[],[],LB,UB,[],options);
            end
            
            
            [~,posLPP]                           =   min(LPP_rep);
            parametersLPP(k_sim,1:numfreeparams) =   parametersLPP_rep(posLPP(1),:);
            LPP(k_sim)                           =   LPP_rep(posLPP(1),:);
            
            con_ll(k_sim, nmodel, con_nmodel)                                 = ll(k_sim);
            con_LPP(k_sim, nmodel, con_nmodel)                                = LPP(k_sim);
            
        end
        
    end
    
    figure;
    set(gcf,'Color',[1,1,1])
    
    
    
    for k = 1:nc
        CC = [k/nc,0,(nc-k)/nc];
        leg{k} = num2str(Ra(k));
        
        subplot(2,2,1)
        hold on
        errorbar(squeeze(mean(A_mat(:,:,k),1)),squeeze(std(A_mat(:,:,k),0,1))./sqrt(n_sims),'-o',...
            'Color',CC,...
            'MarkerFaceColor',CC,...
            'MarkerEdgeColor',CC)
        
        subplot(2,2,2)
        hold on
        errorbar(squeeze(mean(PE_mat(:,:,k),1)),squeeze(std(PE_mat(:,:,k),0,1))./sqrt(n_sims),'-o',...
            'Color',CC,...
            'MarkerFaceColor',CC,...
            'MarkerEdgeColor',CC)
        
        subplot(2,2,3)
        hold on
        plot(suboffers,logitp([Ra(k),Rb(k)],suboffers),':',...
            'LineWidth',1,...
            'Color',CC)
        plot(suboffers,logitp([squeeze(mean(A_mat(:,end,k),1)),mean(Pb_rnd)],suboffers),'-',...
            'LineWidth',2,...
            'Color',CC)
        plot(suboffers,mean(PAsub_mat),'--',...
            'Color',[0,0,0])
        set(gca,'Xlim',[0 10])
        
        subplot(2,2,4)
        hold on
        plot(suboffers,logitp([squeeze(mean(A_mat(:,end,k),1)),mean(Pb_rnd)],suboffers).*(subendow - suboffers),'-',...
            'LineWidth',2,...
            'Color',CC)
        set(gca,'Xlim',[0 10])
        
    end
    
    subplot(2,2,1)
    legend(leg)
    subplot(2,2,2)
    legend(leg)
    subplot(2,2,3)
    errorbar(suboffers,mean(PAsub_mat),std(PAsub_mat)./sqrt(n_sims),'-',...
        'Color',[0,1,0],...
        'MarkerFaceColor',[0,1,0],...
        'MarkerEdgeColor',[0,1,0])
    set(gca,'XLim',[0 10])
    
    
    
    PP = [Px_rnd,Plr1_rnd,Plr2_rnd];           %  Proposer  learning rate
    legB = {'rating temperature','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    for k = 1:3
        
        subplot(2,3,k)
        plot(PP(:,k),parameters(:,k),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(PP(:,k),parameters(:,k));
        
        subplot(2,3,3+k)
        plot(PP(:,k) ,parametersLPP(:,k),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(PP(:,k),parametersLPP(:,k));
        
    end
    
    
end


for k_true = 1:4
    nfpm=[2 3 2 3]; % number of free parameters
    for k_est=1:4
        [~, BIC(:,k_true, k_est)] = aicbic(-con_LPP(:, k_true, k_est), nfpm(k_est), n_trial);
        
        BIC_mean(k_true, k_est)          = mean(BIC(:, k_true, k_est));
        ll_mean(k_true, k_est)           = mean(-con_LPP(:, k_true, k_est));
    end
end
figure
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
title(sprintf('BIC %s', nmodel))
colorbar
