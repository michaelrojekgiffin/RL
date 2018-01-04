clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims  = 100;                           % nsubs to simulates
n_trial = 24;                           % ntrial per cond per session
n_sess  = 2;                            % nsession
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));% parameters of the simulation

% set up conditions and mutliple sessions
%------------------------------------------
Ra      = repmat(-[10,5,0],1,n_sess);            % predator true accepance thereshold (logit intercept)
Rb      = repmat([2,2,2],1,n_sess);              % predator true accpetance moise (logit slope)
n_cond  = size(Ra,2);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
Pa_rnd          = -10 + 10*rand(n_sims,1);%  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 1+2*rand(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
Px_rnd          = 5*rand(n_sims,1);       %  Proposer  rating temperature
Plr_rnd         = 3*rand(n_sims,1);       %  Proposer  learning rate

% setup estimation
%---------------------
options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);
parameters  = NaN(n_sims,2);        parametersLPP  = NaN(n_sims,2);
ll          = NaN(n_sims,1);        LPP            = NaN(n_sims,1);

% Sim loop
%----------
for k_sim = 1:n_sims
    
        % pre-allocate
        O_mat = NaN(n_trial,n_cond);
        D_mat = NaN(n_trial,n_cond);
        
        % get params
        a0 = Pa_rnd(k_sim);
        b0 = Pb_rnd(k_sim);
        bX = Px_rnd(k_sim);
        lr = Plr_rnd(k_sim);
        
    for k_cond = 1:n_cond
        
        % pre-allocat
        %--------------------------
        PA      = NaN(n_trial,numel(offers)); % estimated probability of accepting all offers
        EV      = NaN(n_trial,numel(offers)); % estimated expected value of all offers
        O       = NaN(n_trial,1);             % selected offer, in the euro offer spavce
        Pd      = NaN(n_trial,1);             % trial  probability of accepting the selected offer
        D       = NaN(n_trial,1);             % final decision opf accepting the selected offer
        R       = NaN(n_trial,1);             % Reward
        EG      = NaN(n_trial,1);             % Expected gain
        Pc      = NaN(n_trial,1);             % Estimated prediction of acceptance
        CPE     = NaN(n_trial,1);             % Choice Prediction error
        
        a_t     = NaN(n_trial,1);             % updated logit threshold
        
        % initialize
        %--------------------------
        a_t(1) = a0;
        
        for t = 1:n_trial
            
            % Proposer estimate the decision situation
            %-----------------------------------------------
            PA(t,:)     = logitp([a_t(t,:),b0],offers);     % compute proba of accepting the offers given current model
            EV(t,:)     = (endow - offers).* PA(t,:);   % compute EV of the offers given current model
            
            % Proposer select an Offer
            %-----------------------------------------------
            % Soft max choice (using multinomial choice function)
            p       = exp(bX.*EV(t,:)) ./ sum(exp(bX.*EV(t,:)));        % multinomial choice function
            pd      = makedist('multinomial','probabilities',p);        % estimate the pdf from the pC
            ypdf    = pdf(pd,1:numel(offers));                          % generate pdf for the offers
            kO      = random(pd);                                       % selected offer, in the 1:numel(offer) spavce
            O(t)    = offers(kO);                                       % resample Offer in pdf (="soft-max")
            
            % Proposer make choices and observe decision
            %-------------------------------------------
            Pd(t)   = logitp([Ra(k_cond),Rb(k_cond)],O(t)); % Reciever Estimated accepantce proba of accepting the offer
            D(t)    = double(rand(1)<Pd(t));                % Sampling Reciever's decision given the proba.
            Pc(t)   = logitp([a_t(t,:),b0],O(t));
            
            % Updating Proposer estimation of the reciever's acceptance
            % function
            %------------------------------------------------------------
            CPE(t)      = D(t) - Pc(t);
            a_t(t+1,:)  = a_t(t,:) + lr.*CPE((t)); % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
            
        end
        
        O_mat(:,k_cond) = O;
        D_mat(:,k_cond) = D;
        
    end
    
    
    n_rep           = 10;
    parameters_rep  = NaN(n_rep,2);     parametersLPP_rep  = NaN(n_rep,2);
    ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
    
    lb = [-15 0];       LB = [-Inf 0];
    ub = [5 5];         UB = [Inf Inf];
    ddb = ub - lb;
    
    for k_rep = 1:n_rep
        x0 = lb + rand(1,2).*ddb;
        %standard estimation
        [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) learning_model_2params(x,O,D,a0,b0),x0,[],[],[],[],lb,ub,[],options);
        %lalace estimation
        [parametersLPP_rep(k_rep,1:2),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,a0,b0),x0,[],[],[],[],LB,UB,[],options);
    end
    
    [~,pos] = min(ll_rep); 
    parameters(k_sim,:)    =   parameters_rep(pos(1),:);
    ll(k_sim)              =   ll_rep(pos(1),:);
    
    
    [~,posLPP] = min(LPP_rep); 
    parametersLPP(k_sim,:)    =   parametersLPP_rep(posLPP(1),:);
    LPP(k_sim)              =   LPP_rep(posLPP(1),:);
    
end


figure;
set(gcf,'Color',[1,1,1])
subplot(2,2,1)
plot(Px_rnd ,parameters(:,1),'o',...
    'MarkerEdgeColor',[0,0,0],...
    'MarkerFaceColor',[1,1,1])
xlabel('true rating temperature')
ylabel('estimated rating temperature')

[corrR(1),corrP(1)] = corr(Px_rnd ,parameters(:,1));

subplot(2,2,2)
plot(Plr_rnd ,parameters(:,2),'o',...
    'MarkerEdgeColor',[0,0,0],...
    'MarkerFaceColor',[1,1,1])
xlabel('true learning rate')
ylabel('estimated learning rate')

[corrR(2),corrP(2)] = corr(Plr_rnd ,parameters(:,2));



subplot(2,2,3)
plot(Px_rnd ,parametersLPP(:,1),'o',...
    'MarkerEdgeColor',[0,0,0],...
    'MarkerFaceColor',[1,1,1])
xlabel('true rating temperature')
ylabel('estimated rating temperature LPP')

[corrR_LPP(1),corrP_LPP(1)] = corr(Px_rnd ,parametersLPP(:,1));

subplot(2,2,4)
plot(Plr_rnd ,parametersLPP(:,2),'o',...
    'MarkerEdgeColor',[0,0,0],...
    'MarkerFaceColor',[1,1,1])
xlabel('true learning rate')
ylabel('estimated learning rate LPP')

[corrR_LPP(2),corrP_LPP(2)] = corr(Plr_rnd ,parametersLPP(:,2));