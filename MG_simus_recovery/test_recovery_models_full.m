clear
close all
clc

resp_params = load('resp_params');
resp_params = resp_params.resp_params([1, 3, 5, 2, 4], :);
% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_trial_per_block = 120;                % the number of trials per block that all conditions will be within = 120

n_sims  = 40;                           % nsubs to simulates
% n_trial = 12;                           % ntrial per cond per session
% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
%     Pa_rnd          = -9 + 6*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
%     Pb_rnd          = 2.9+.2*rand(n_sims,1);  %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
%
Pa_rnd          = -3 + 2*rand(n_sims,1);  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = .2+.5*rand(n_sims,1);  %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)

Px_rnd          = .5+2.5*rand(n_sims,1);  %  Proposer  rating temperature
% Px_rnd          = 3+3*rand(n_sims,1);   %  Proposer  rating temperature
Plr1_rnd        = rand(n_sims,1);         %  Proposer  learning rate
Plr2_rnd        = rand(n_sims,1);         %  Proposer  learning rate

for condsim = 3:5 % amount of conditions to include
    % n_sess determines how many conditions the simulation plays against,
    % it's determined in Ra and Rb, both of whose dimensions are used to
    % determine n_cond in learning_models_timeseries and learning_models_estim
    n_sess  = condsim;                      % number of different opponents I'm playing against
    
    n_trial = round(n_trial_per_block/condsim);    % ntrial per cond per session
    
    %     cond2learn  = -[12,9,6,3,0];
    cond2learn  = resp_params(1:condsim, 1)';
    
    offers  = 0:20;
    endow   = 20*ones(1,numel(offers));% parameters of the simulation
    
    modelspace = [1 2 3 4];
    nfpm=[2 3 2 3];
    nmods = numel(modelspace);
    
    % set up conditions and mutliple sessions
    %------------------------------------------
    % cond2learn  = -[12,9,6,3,0];
    nc          = numel(cond2learn);
    Ra          = repmat(cond2learn,1,n_sess);           % responder true accepance thereshold (logit intercept)
    %     Rb          = repmat(3*ones(1,nc),1,n_sess);         % responder true acceptance noise (logit slope)
    Rb          = repmat(.4*ones(1,nc),1,n_sess);         % responder true acceptance noise (logit slope)
    n_cond      = size(Ra,2);
    
    
    % setup estimation
    %---------------------
    options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000, 'display', 'off');
    parameters  = NaN(n_sims,3,nmods,nmods);        parametersLPP  = NaN(n_sims,3,nmods,nmods);
    ll          = NaN(n_sims,nmods,nmods);        LPP            = NaN(n_sims,nmods,nmods);
    
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
            
            [O,D] = learning_models_timeseries([bX,lr1,lr2],[Ra;Rb],n_trial,a0,b0,ktm);
            
            
            lb = [0 0 0];          LB = [0 0 0];
            ub = [15 1 1];         UB = [Inf 10 1];
            ddb = ub - lb;
            
            for kem = modelspace
                fprintf('running true model %d, estimated model %d, %d conditions, simulation %d out of %d\n', ktm, kem, condsim, k_sim, n_sims);
                
                n_rep           = 5;
                parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
                ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
                
                for k_rep = 1:n_rep
                    x0 = lb + rand(1,3).*ddb;
                    % %standard estimation
                    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,a0,b0,kem),x0,[],[],[],[],lb,ub,[],options);
                    % %lalace estimation
                    [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2(x,O,D,a0,b0,kem),x0,[],[],[],[],LB,UB,[],options);
                end
                [~,pos] = min(ll_rep);
                parameters(k_sim,:,ktm,kem)    =   parameters_rep(pos(1),:);
                ll(k_sim,ktm,kem)              =   ll_rep(pos(1),:);
                
                [~,posLPP] = min(LPP_rep);
                parametersLPP(k_sim,:,ktm,kem)      =   parametersLPP_rep(posLPP(1),:);
                LPP(k_sim,ktm,kem)                  =   LPP_rep(posLPP(1),:);
            end
        end
    end
    
    time = clock;
    time = strcat(num2str(time(4)), num2str(time(5)));
    save(['MG_recovery_', num2str(condsim),'_', date])
end

fprintf('FINISHED!!!!\n');
%%
for k_true = modelspace
    
    MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
    LL = squeeze(ll(:,k_true,:));
    
    for k_est= modelspace
        bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
    end
    [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
    % [postBMC,outBMC]=VBA_groupBMC(-LL');
    BMC_output(k_true).post = postBMC;
    BMC_output(k_true).out = outBMC;
    
    Ep(k_true,:) = 100*BMC_output(k_true).out.ep;
    
end


figure
set(gcf,'Color',[1,1,1])


colormap(flipud(gray))
imagesc(flipud(Ep))
ylabel('simulated model #')
xlabel('estimated model #')
set(gca,'XTick',1:4,...
    'YTick',1:4,...
    'XTickLabel',(1:4),...
    'YTickLabel',fliplr(1:4))

c = colorbar;
c.Label.String = 'Exceedance probability (%)';


for k_true = modelspace
    
    
    legB = {'rating temperature','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    
    
    title(strcat(['Model ',num2str(k_true)]));
    for k = 1:3
        
        subplot(2,3,k)
        plot(MP(:,k),squeeze(parameters(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(MP(:,k),squeeze(parameters(:,k,k_true,k_true)));
        
        subplot(2,3,3+k)
        plot(MP(:,k) ,squeeze(parametersLPP(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),squeeze(parametersLPP(:,k,k_true,k_true)));
        
    end
    
end