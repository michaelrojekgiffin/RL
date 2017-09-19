% clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:1:10;
ntrial  = 30;
endow   = 10*ones(1,numel(offers));

% parameters of the simulation
%--------------------------
v0      = -5;       %  Prevy initial prior on thereshold
B0      = 3;       %  Prey rating temperature
a0      = 1;       %  Prey learning rate 
B       = 1;       %  Prey learning rate 

v1      = -15;        % predator true thereshold
B1      = 2;        % predator true thereshold


% logistic choice function
%--------------------------
% logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));

logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));

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




figure;
hold on
plot(V,'r');
plot(O,'b');
plot(D,'k')

legend('V','O','D')


figure;
subplot(1,2,1)
hold on
plot(offers,logitp([v1,B1],offers),'r');
plot(offers,logitp([V(1),B],offers),'b');
plot(offers,logitp([V(end),B],offers),'g');
legend('true','prior','estimated')
title('acceptance function')


subplot(1,2,2)
hold on
plot(offers,EV(1,:),'b');
plot(offers,EV(end,:),'g');
legend('prior','estimated')
title('EV function')