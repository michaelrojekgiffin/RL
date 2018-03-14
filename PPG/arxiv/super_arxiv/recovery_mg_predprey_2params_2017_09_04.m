clear
close all
% clc

predprey = 'predator';

switch predprey
    case 'predator'
        priors = load('data_matlab/predator_priors.mat');
    case 'prey'
        priors = load('data_matlab/prey_priors.mat');
end

% randomize generator seed
%--------------------------
rng('shuffle')

nsims = 20;
% parameters of the task
%--------------------------
offers  = 0:1:10;
ntrial  = 10;
endow   = 10*ones(1,numel(offers));

% parameters of the simulation
%--------------------------

v1      = -15;      % receiver (opponent) true thereshold, intercept
B1      = 2.5;        % reciever (opponent) true thereshold, slope


% logistic choice function
%--------------------------
% logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));

logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));

v0 = priors.parameters(1); %  proposer initial prior on threshold (estimated from concatenated first trial investments) 
B = priors.parameters(2);   %  proposer noise (proposer's estimate of slope of opponent, (estimated from concatenated first trial investments) 

B0rnd      = 10*rand(nsims,1);       %  proposer rating temperature
a0rnd      = 10*rand(nsims,1);       %  proposer learning rate

options     = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);


for k_sim = 1:nsims
    
    B0 = B0rnd(k_sim);  % proposer temperature 
    a0 = a0rnd(k_sim);  % proposer learning rate

    % pre-allocat
    %--------------------------
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
    V       = NaN(ntrial,1);             % selected offer
    kO      = NaN(ntrial,1);             % selected offer, in the 1:numel(offer) spavce
    O       = NaN(ntrial,1);             % selected offer, in the euro offer spavce
    Pd      = NaN(ntrial,1);             % trial  probability of accepting the selected offer
    D       = NaN(ntrial,1);             % final decision opf accepting the selected offer
    R       = NaN(ntrial,1);             % Reward
    EG      = NaN(ntrial,1);             % Expected gain
    CPE     = NaN(ntrial,1);             % Choice Prediction error
    
    
    % initialize
    %--------------------------
    V(1) = v0;
    
    for t = 1:ntrial
        
        % Proposer estimate the decision situation
        %-----------------------------------------------
        % PA(t,:)     = maxp(V(t,:),offers);         % compute proba of accepting the offers given current model
        PA(t,:)     = logitp([V(t,:),B],offers);     % compute proba of accepting the offers given current model
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
        p = exp(B0.*EV(t,:)) ./ sum(exp(B0.*EV(t,:)));          % multinomial choice function
        pd = makedist('multinomial','probabilities',p);         % estimate the pdf from the pC
        ypdf = pdf(pd,1:numel(offers));                         % generate pdf for the offers
        kO(t) = random(pd);                                     % select index of pd based on probability created by makedist
        O(t) = offers(kO(t));                                   % resample Offer in pdf (="soft-max")
        
        % Proposer make choices and observe decision
        %-------------------------------------------
        Pd(t)       = logitp([v1,B1],O(t));                     % Reciever accepantce proba of accepting the offer, given true 
                                                                % values of receivers acceptance function
        
        % D is whether the investment (proposal of our simulated subject)
        % is "accepted" (kill/survive) or "rejected" (miss/die)
        D(t)        = double(rand(1)<Pd(t));              % Sampling Reciever's decision given the proba.
        
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
    
    n_rep           = 5;
    parameters_rep  = NaN(n_rep,2);
    ll_rep          = NaN(n_rep,1);
    
    lb = [0 0];
    ub = [50 10];
    ddb = ub - lb;
    
    for k_rep = 1:n_rep     
        x0 = lb + rand(1,2).*ddb;
        [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_2params_2017_09_04(x,[v0, B],O,D,1,predprey),x0,[],[],[],[],lb,ub,[],options);
    end
    
    [~,pos] = min(ll_rep);
    % parameters(k_sim,:)    =   parameters_rep;
    parameters(k_sim,:)    =   parameters_rep(pos(1),:);
    ll(k_sim)              =   ll_rep(pos(1),:);
    
end

figure
plot(offers,logitp([v1,B1],offers),'r');

figure;
subplot(2,1,1)
plot(B0rnd,parameters(:,1),'o')
[RR, PP] = corrcoef(B0rnd,parameters(:,1));
txt1 = sprintf('r = %f\n p = %f', RR(2), PP(2));
text(mean(B0rnd), mean(parameters(:,1)), txt1);
title('proposer temperature');
xlabel('true parameters')
ylabel('estimated parameters')
lsline;

subplot(2,1,2)
plot(a0rnd,parameters(:,2),'o')
[RR, PP] = corrcoef(a0rnd,parameters(:,2));
txt1 = sprintf('r = %f\n p = %f', RR(2), PP(2));
text(mean(a0rnd), mean(parameters(:,2)), txt1);
title('proposer alpha (learning rate)');
xlabel('true parameters')
ylabel('estimated parameters')
lsline;
% 
% subplot(2,2,2)
% plot(v0rnd,parameters(:,2),'o')
% [RR, PP] = corrcoef(v0rnd,parameters(:,2));
% txt1 = sprintf('r = %f\n p = %f', RR(2), PP(2));
% text(mean(v0rnd), mean(parameters(:,2)), txt1);
% title('proposer estimate of receiver intercept');
% xlabel('true parameters')
% ylabel('estimated parameters')
% lsline;

% 
% subplot(2,2,4)
% plot(Brnd,parameters(:,4),'o')
% [RR, PP] = corrcoef(Brnd,parameters(:,4));
% txt1 = sprintf('r = %f\n p = %f', RR(2), PP(2));
% text(mean(Brnd), mean(parameters(:,4)), txt1);
% title('proposer estimate of receiver slope');
% xlabel('true parameters')
% ylabel('estimated parameters')
% lsline;