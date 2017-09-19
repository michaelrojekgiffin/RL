clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:.1:10;
ntrial  = 12;
endow   = 10*ones(1,numel(offers));

% parameters of the simulation
%--------------------------
bR      = [5,2];        % true receiver B
bhat0   = [0,.5];       % Proposer prior on receiver B (mean)
shat0   = [3,3];        % Proposer prior on receiver B (std)

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));

% pre-allocate
%--------------------------
PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
pos     = NaN(ntrial,1);             % estimated probability of accepting all offers
V       = NaN(ntrial,1);             % selected offer
Vmax    = NaN(ntrial,1);             % selected offer, calculated with greedy max
pV      = NaN(ntrial,1);             % selected offer, calculated with soft max
Pr      = NaN(ntrial,1);             % trial estimated probability of accepting the selected offer
D       = NaN(ntrial,1);             % observed decision
MSE     = NaN(ntrial,1);             % Mean Square Error
B = .5;                              % multinomial choice function temperature

% initialize
%--------------------------
bhat = bhat0;
shat = shat0;

for t = 1:ntrial
    
    % check the distance between true and estimated accepatnce function
    %--------------------------------------------------------------------
    MSE(t) = sum((logitp(bhat(t,:),offers) - logitp(bR,offers)).^2)./numel(offers);
    
   % Proposer estimate the decision situation
   %-----------------------------------------------
    PA(t,:)     = logitp(bhat(t,:),offers);     % compute proba of accepting the offers given current model
    EV(t,:)     = (endow - offers).* PA(t,:);   % compute EV of the offers given current model
    
   % Proposer select an Offer
   %-----------------------------------------------
    % Option1: Hard max choice
    [~,pos(t)]  = max(EV(t,:));
    Vmax(t)     = offers(pos(t));
    
    % Option2: Soft max choice (using multinomial choice function)
    p = exp(B.*EV(t,:)) ./ sum(exp(B.*EV(t,:)));            % multinomial choice function
    pd = makedist('multinomial','probabilities',p);         % estimate the pdf from the pC
    ypdf = pdf(pd,1:numel(offers));                         % generate pdf for the offers
    pV(t) = offers(random(pd));                             % resample choice in pdf (="soft-max")
    
    %     figure;hold on; plot(offers,p,'ob');plot(offers,ypdf,'r') % plot
    
    % Proposer make choices and observe decision
    %-------------------------------------------
    V(t)        = pV(t);                                    % Select Offer. Choose here pV(t) or Vmax(t)
    Pr(t)       = logitp(bR,V(t));                          % Reciever true accepantce proba of accepting the offer
    D(t)        = double(rand(1)<Pr(t));                    % Sampling Reciever's decision given the proba.
    
    % Updating Proposer estimation of the reciever's acceptance
    % function
    %------------------------------------------------------------
    prior1      = @(b1) normpdf(b1,bhat(t,1),shat(t,1));    % prior for intercept
    prior2      = @(b2) normpdf(b2,bhat(t,2),shat(t,2));    % prior for slope
    
    % combine prior and likelihood
    post        = @(b) prod(binopdf(D(t),ones(1),logitp(b,V(t)))) ...  % likelihood
        * prior1(b(1)) * prior2(b(2));   
    
    % MCM estimation of the posterior
    initial = bhat0;
    nsamples = 500;
    trace = slicesample(initial,nsamples,'pdf',post,'burnin',100);
    
    % figure
    % subplot(2,1,1);plot(trace(:,1));ylabel('Intercept');xlabel('Sample Number');
    % subplot(2,1,2);plot(trace(:,2));ylabel('Slope');xlabel('Sample Number');
    
    % Update: posterior becomes prior
    bhat(t+1,:) = mean(trace);
    shat(t+1,:) = std(trace);
   
end

%==========================%
% Figure 1
%==========================%
h1 = figure('Units', 'pixels', ...
    'Position', [400 300 800 500]);
set(gcf,'Color',[1,1,1])

% plot 1: compare true, initial estimate and final estimate of R-accepatnce function 
subplot(2,2,1)
hold on
plot(offers,logitp(bR,offers),'-r')
plot(offers,logitp(bhat0,offers),'-.g')
plot(offers,logitp(bhat(end,:),offers),'--b')
legend('true \phi','prior \phi','posterior \phi')
xlabel('Offer (euro)')
ylabel('p(Accept)')

% Plot 1: "Learning": difference between estimated and true acceptance
% function
subplot(2,2,2)
hold on
plot(MSE,'-ok',...
    'MarkerFaceColor',[0,0,0])
ylabel('MSE')
xlabel('trial')
legend('distance true vs estimated \phi')

% plot 3: compare true, initial estimate and final estimate of offers' expected values 
subplot(2,2,3)
hold on
[Opt,pos] = max((endow - offers).*logitp(bR,offers));
Voptimal = offers(pos);
[OptHat,poshat] = max((endow - offers).*logitp(bhat(end,:),offers));
Vhat = offers(poshat);
[Opt0,pos0] = max((endow - offers).*logitp(bhat0,offers));
V0 = offers(pos0);

plot(offers,(endow - offers).*logitp(bR,offers),'-r')
plot(offers,(endow - offers).*logitp(bhat0,offers),'-.g')
plot(offers,(endow - offers).*logitp(bhat(end,:),offers),'--b')

plot([Voptimal Voptimal],[0 Opt],'-r')
plot([V0 V0],[0 Opt0],'-.g')
plot([Vhat Vhat],[0 OptHat],'--b')
legend('actual','prior (initial)','estimated (final)')
xlabel('Offer (euro)')
ylabel('Expected value')

% plot 3: Display chosen values (selected with greedy max of soft max)
subplot(2,2,4)
hold on
plot(Vmax,'-ok',...
    'MarkerFaceColor',[0,0,0])
plot(pV,'-ob',...
    'MarkerFaceColor',[0,0,1])
plot([1 ntrial],[Voptimal Voptimal],'--r')
set(gca,'YLim',[min(offers) max(offers)])
xlabel('trial')
ylabel('euro')
legend('Offer Max','Offer Stoch','Optimal Offer')


%==========================%
% Figure 2
%==========================%
h2 = figure('Units', 'pixels', ...
    'Position', [300 400 1000 300]);

set(gcf,'Color',[1,1,1])

% plot 1: compare true, initial estimate and final estimate parameters distributions 
subplot(1,3,1)
xx = -10:.1:10;
D1 = normpdf(xx,repmat(bhat0(1),1,numel(xx)),repmat(shat0(1),1,numel(xx)));
D2 = normpdf(xx,repmat(bhat0(2),1,numel(xx)),repmat(shat0(2),1,numel(xx)));
Da = normpdf(xx,repmat(bhat(end,1),1,numel(xx)),repmat(shat(end,1),1,numel(xx)));
Db = normpdf(xx,repmat(bhat(end,2),1,numel(xx)),repmat(shat(end,2),1,numel(xx)));
hold on
plot([bR(1),bR(1)],[0 .5],'-','Color',[.5,.5,0],'LineWidth',2)
plot([bR(2),bR(2)],[0 .5],'-','Color',[0,.5,.5],'LineWidth',2)
plot(xx,D1,':','Color',[.5,.5,0])
plot(xx,D2,':','Color',[0,.5,.5])
plot(xx,Da,'--','Color',[.5,.5,0])
plot(xx,Db,'--','Color',[0,.5,.5])
legend('true \beta_0','true \beta_1',...
    'prior \beta_0','prior \beta_1',...
    'posterior \beta_0','posterior \beta_1')

% Plot 2: Evolution of parameter means
subplot(1,3,2)
hold on
plot(bhat(:,1),'-ob',...
    'MarkerFaceColor',[0,0,1])
plot(bhat(:,2),'-or',...
    'MarkerFaceColor',[1,0,0])
plot([1 ntrial],[bR(1) bR(1)],'--b')
plot([1 ntrial],[bR(2) bR(2)],'--r')
legend('\beta_0','\beta_1')
ylabel('Parameter mean')
xlabel('trial')

% Plot 3: Evolution of parameter std
subplot(1,3,3)
hold on
plot(shat(:,1),'-ob',...
    'MarkerFaceColor',[0,0,1])
plot(shat(:,2),'-or',...
    'MarkerFaceColor',[1,0,0])
legend('\beta_0','\beta_1')
ylabel('Parameter std')
xlabel('trial')
