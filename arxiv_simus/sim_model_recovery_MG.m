% clear
close all
% clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:1:10;
ntrial  = 60;
endow   = 10*ones(1,numel(offers));

options     = optimset('Algorithm', 'interior-point', 'Display', 'off', 'MaxIter', 10000);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));

% pre-allocat
%--------------------------
PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
V       = NaN(ntrial,1);             % selected offer
kO      = NaN(ntrial,1);            % selected offer, in the 1:numel(offer) spavce
O       = NaN(ntrial,1);             % selected offer, in the euro offer spavce
Pd      = NaN(ntrial,1);             % trial  probability of accepting the selected offer
D       = NaN(ntrial,1);             % final decision of accepting the selected offer
R       = NaN(ntrial,1);             % Reward
EG      = NaN(ntrial,1);             % Expected gain
PE      = NaN(ntrial,1);             % Choice Prediction error


% parameters of the simulation
%--------------------------
v0_array      = 0:10;       %  Prey initial prior on thereshold
B0_array      = 0:.5:10;       %  Prey rating temperature - range
a0_array      = 0.05:0.05:1;       %  Prey learning rate - range


% to examine recovery
%--------------------------
N_sim = 30;             % how many simulations to run
sim_rec = NaN(N_sim, 6); % matrix representing recovered (cols 1:3) and true parameters (cols 4:6)

% pre-allocate
%--------------------------

v1      = 8;        % predator true thereshold
B1      = .01;        % predator true thereshold

% initialize
%--------------------------
for ss = 1:N_sim
    % run simiulation with random combinations of parameteres
    V(1) = datasample(v0_array, 1);
    a0 = datasample(a0_array, 1);
    B0 = datasample(B0_array, 1);
    
    for t = 1:ntrial
        % Proposer estimate the decision situation
        %-----------------------------------------------
        PA(t,:)     = maxp(V(t,:),offers);     % compute proba of accepting the offers given current model
        EV(t,:)     = (endow - offers).* PA(t,:);   % compute EV of the offers given current model
        
        % Proposer select an Offer
        %-----------------------------------------------
        % Soft max choice (using multinomial choice function)
        p = exp(B0.*EV(t,:)) ./ sum(exp(B0.*EV(t,:)));          % multinomial choice function
        pd = makedist('multinomial','probabilities',p);         % estimate the pdf from the pC
        ypdf = pdf(pd,1:numel(offers));                         % generate pdf for the offers
        kO(t) = random(pd);
        O(t) = offers(kO(t));                                   % resample Offer in pdf (="soft-max")
        
        % Proposer make choices and observe decision
        %-------------------------------------------
        Pd(t)       = logitp([v1/2,B1],O(t));               % Reciever Estimated accepantce proba of accepting the offer
        D(t)        = double(rand(1)<Pd(t));                % Sampling Reciever's decision given the proba.
        
        if D(t) == 1
            R(t) = 10 - O(t);
        else
            R(t) = 0;
        end
        % Updating Proposer estimation of the reciever's acceptance
        % function
        %------------------------------------------------------------
        %                 PE(t) = EG(t) - R(t);
        
        V(t+1,:) = V(t,:) - a0.*sign(D(t)-.5); % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        
    end

    sub_o               = O;
    sub_r               = R;
    
    n_rep               = 5;
    parameters_rep      = NaN(n_rep,3);
    ll_rep              = NaN(n_rep,1);
    
    % IMPORTANT - I changed the upper bound on the third parameter (alpha)
    % from 10 to 1
    for k_rep = 1:n_rep
        [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_Full_2016_11(x,sub_o,sub_r,1),[10*rand() 10*rand() rand() ],[],[],[],[],[0 0 0],[10 10 1],[],options);
    end
    
    [~,pos]             = min(ll_rep);
    
    parameters          =   parameters_rep(pos(1),:); % parameters are B0, V0, and a0, respectively
    %             ll(rec_count)              =   ll_rep(pos(1),:);
    
    sim_rec(ss, 1:3)    = parameters;                 % store recovered parameters
    sim_rec(ss, 4:6)    = [B0, V(1), a0];             % store true parameters
end


figure;
scatter(sim_rec(:, 1), sim_rec(:, 4)); 
lsline;
title('Recovered Beta vs Real beta');
xlabel('Recovered');
ylabel('Real');

figure;
scatter(sim_rec(:, 2), sim_rec(:, 5)); 
lsline;
title('Recovered V0 vs Real V0');
xlabel('Recovered');
ylabel('Real');


figure;
scatter(sim_rec(:, 3), sim_rec(:, 6)); 
lsline;
title('Recovered alpha vs Real alpha');
xlabel('Recovered');
ylabel('Real');
