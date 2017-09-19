clear
close all
% clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
offers  = 0:.1:10;
ntrial  = 30;
endow   = 10*ones(1,numel(offers));

% parameters of the simulation
%--------------------------
% vR      = 3;       % true receiver B1
% vS      = 1;       % true receiver b2
% VP      = -2;       % Proposer initial prior on receiver B (mean)
% VS      = 4;       % Proposer temperature receiver B (mean)

vR      = -6;       % true receiver B1
vS      = 2;       % true receiver b2
VP      = 3;       % Proposer initial prior on receiver B (mean)
VS      = 1;       % Proposer temperature receiver B (mean)


a0      = 1;
a1      = 3;

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x-median(offers)))./(1+exp(b(1)+b(2).*(x-median(offers))));
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));

% pre-allocat
%--------------------------
PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
pos     = NaN(ntrial,1);             % estimated probability of accepting all offers
V       = NaN(ntrial,1);             % selected offer
Vmax    = NaN(ntrial,1);             % selected offer, calculated with greedy max
pV      = NaN(ntrial,1);             % selected offer, calculated with soft max
Pd      = NaN(ntrial,1);             % trial estimated probability of accepting the selected offer
D       = NaN(ntrial,1);             % observed decision
CPE     = NaN(ntrial,1);             % Choice Prediction error
trPE    = NaN(ntrial,1);             % Transitivity Prediction Error
MSE     = NaN(ntrial,1);             % Mean Square Error
B = .5;                              % multinomial choice function temperature

% initialize
%--------------------------
bR      = [vR,vS];
bhat0   = [VP,VS];
bhat    = bhat0;

for t = 1:ntrial
    
    % check the distance between true and estimated accepatnce function
    %--------------------------------------------------------------------
    % MSE(t) = sum((logitp(bhat(t,:),offers) - maxp(VR,offers)).^2)./numel(offers);
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
    Pd(t)       = logitp(bhat(t,:),V(t));                   % Reciever Estimated accepantce proba of accepting the offer
    D(t)        = double(rand(1)<logitp(bR,V(t)) );         % Sampling Reciever's decision given the proba.
    % D(t)        = maxp(VR,V(t));                           % Sampling Reciever's decision .
    L(t) = prod(binopdf(D(t),ones(1),logitp(bhat(t,:),V(t))));
    
    % Updating Proposer estimation of the reciever's acceptance
    % function
    %------------------------------------------------------------
    CPE(t) = D(t) - Pd(t);
    bhat(t+1,1) = bhat(t,1) + a0*CPE(t);
    
    
    nrep = 100;
    if t <2
        trPE(t) = 0;
    else
        C = nchoosek(1:t,2);
        ntrI = mean(V(C(:,1))>=V(C(:,2)) == (D(C(:,1))>=D(C(:,2))));            % count the number of "intransitive choices" in history
        
        randC = rand(t*nrep,1);                                                 % setup random Nature
        modpC = repmat(logitp([bhat(t+1,1),bhat(t,2)],V(1:t)),nrep,1);          % simulate history of decisions probabilitiees given current params and history of offers
        Dsims = double(randC<modpC);                                            % compute simulated decisions given randomNature
        repC = repmat(C,nrep,1);
        ntrIsims= mean(V(repC(:,1))>=V(repC(:,2)) == (Dsims(repC(:,1))>=Dsims(repC(:,2)))); % count the number of "intransitive choices" in this sim data

        trPE(t) = ntrI-ntrIsims;
    end
    
    bhat(t+1,2) = bhat(t,2) + a1.*trPE(t);
    
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
% plot(offers,maxp(bR,offers),'-r')
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
%[Opt,pos] = max((endow - offers).*maxp(bR,offers));
[Opt,pos] = max((endow - offers).*logitp(bR,offers));
Voptimal = offers(pos);
[OptHat,poshat] = max((endow - offers).*logitp(bhat(end,:),offers));
Vhat = offers(poshat);
[Opt0,pos0] = max((endow - offers).*logitp(bhat0,offers));
V0 = offers(pos0);

% plot(offers,(endow - offers).*maxp(bR,offers),'-r')
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

% Plot 2: Evolution of parameter means
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

