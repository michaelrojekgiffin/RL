 clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

nsims = 20;
% parameters of the task
%--------------------------
offers  = 0:1:10;
ntrial  = 60;
endow   = 10*ones(1,numel(offers));

% parameters of the simulation
%--------------------------



v1      = -15;        % predator true thereshold
B1      = 2;        % predator true thereshold


% logistic choice function
%--------------------------
% logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));

logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));


v0rnd      = -10 + 10*rand(nsims,1);%  Prey initial prior on threshold (to be fitted or estimated from no-feedback games)
B0rnd      = 5*rand(nsims,1);       %  Prey rating temperature
a0rnd      = 3*rand(nsims,1);       %  Prey learning rate

Brnd       = 5*rand(nsims,1);       %  Prey noise (to be fitted or estimated from no-feedback games)


options     = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);


for k_sim = 1:nsims
    
    v0 = v0rnd(k_sim);
    B0 = B0rnd(k_sim);
    a0 = a0rnd(k_sim);
    B = Brnd(k_sim);
    
    % pre-allocat
    %--------------------------
    PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
    EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
    V       = NaN(ntrial,1);             % selected offer
    kO       = NaN(ntrial,1);            % selected offer, in the 1:numel(offer) spavce
    O       = NaN(ntrial,1);             % selected offer, in the euro offer spavce
    Pd      = NaN(ntrial,1);             % trial  probability of accepting the selected offer
    D      = NaN(ntrial,1);             % final decision opf accepting the selected offer
    R       = NaN(ntrial,1);             % Reward
    EG      = NaN(ntrial,1);             % Expected gain
    CPE      = NaN(ntrial,1);             % Choice Prediction error
    
    
    
    
    % initialize
    %--------------------------
    V(1) = v0;
    
    for t = 1:ntrial
        
        % Proposer estimate the decision situation
        %-----------------------------------------------
        % PA(t,:)     = maxp(V(t,:),offers);     % compute proba of accepting the offers given current model
        PA(t,:)     = logitp([V(t,:),B],offers);     % compute proba of accepting the offers given current model
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
        Pd(t)       = logitp([v1,B1],O(t));               % Reciever Estimated accepantce proba of accepting the offer
        D(t)        = double(rand(1)<Pd(t));                % Sampling Reciever's decision given the proba.
        Pc(t)       = logitp([V(t,:),B],O(t));
        
        
        % Updating Proposer estimation of the reciever's acceptance
        % function
        %------------------------------------------------------------
        CPE(t) = D(t) - Pc(t);
        V(t+1,:) = V(t,:) + a0.*CPE((t)); % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        
    end
%     
%     figure;
%     hold on
%     plot(V,'r');
%     plot(O,'b');
%     plot(D,'k')
%     legend('V','O','D')

    
    n_rep           = 5;
    parameters_rep  = NaN(n_rep,4);
    ll_rep          = NaN(n_rep,1);
    
    lb = [0 -15 0 0];
    ub = [10 5 5 10];
    ddb = ub - lb;
    
    for k_rep = 1:n_rep     
        x0 = lb + rand(1,4).*ddb;
        [parameters_rep(k_rep,1:4),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation3_Full_2016_11(x,O,D,1),x0,[],[],[],[],lb,ub,[],options);
    end
    
    [~,pos] = min(ll_rep);
    % parameters(k_sim,:)    =   parameters_rep;
     parameters(k_sim,:)    =   parameters_rep(pos(1),:);
    ll(k_sim)              =   ll_rep(pos(1),:);
    
end

figure
plot(offers,logitp([v1,B1],offers),'r');

figure;
subplot(2,2,1)
plot(B0rnd,parameters(:,1),'o')
lsline;
subplot(2,2,2)
plot(v0rnd,parameters(:,2),'o')
lsline;
subplot(2,2,3)
plot(a0rnd,parameters(:,3),'o')
lsline;
subplot(2,2,4)
plot(Brnd,parameters(:,4),'o')
lsline;