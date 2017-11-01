clear
close all
clc
 


offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));% parameters of the simulation


%% level-0
% at level 0, players play at random

opp =  makedist('uniform','lower',offers(1),'upper',offers(end));           % estimate the pdf from the pC
X = pdf(opp, offers);
figure;
set(gcf,'Color',[1,1,1])
plot(offers,X,'-o',...
    'MarkerFaceColor',[1,1,1])
title('Level 0')
set(gca,'YLim',[0 .2])
xlabel('offers')
ylabel('likelihood')


%% level-1
% at level 1, players expect opponents to be level-0 player

opp =  makedist('Uniform','lower',offers(1),'upper',offers(end));           % estimate the pdf from the pC
X = pdf(opp, offers);
Y = cumsum(X)-X(1);

EVprey     = (endow - offers).* (Y);                       % compute EV of the offers given current model
% EVpred     = (endow - offers) + (endow - offers).* (Y);   % compute EV of the offers given current model
EVpred = NaN(1, length(offers));
for tc = 1:length(EVpred)
    EVpred(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (X(1:tc)));
end
EVpred(1) = 10; % EV for an investment of 0 should always be 10


PCprey     = EVprey./sum(EVprey);                       % compute EV of the offers given current model
PCpred     = EVpred./sum(EVpred);   

figure
set(gcf,'Color',[1,1,1])

subplot(2,2,1)
plot(offers,X,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('opponent probability')
set(gca,'YLim',[0 .2])
title('Level 1')

subplot(2,2,2)
plot(offers,Y,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('p(X>(opponent))')

subplot(2,2,3)
hold on
plot(offers,EVprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,EVpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('EV')
legend('prey','predator')

subplot(2,2,4)
hold on
plot(offers,PCprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,PCpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('likelihood')
legend('prey','predator')


%% level-2
% at level 2, players expect opponent to be level-1 players

XPrey = PCprey;
XPred = PCpred;

YPrey = cumsum(XPrey);
YPred = cumsum(XPred);

EVprey     = (endow - offers).* (YPred);                       % compute EV of the offers given current model
% EVpred     = (endow - offers) + (endow - offers).* (Y);   % compute EV of the offers given current model
EVpred = NaN(1, length(offers));
for tc = 1:length(EVpred)
    EVpred(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (XPred(1:tc)));
end
EVpred(1) = 10; % EV for an investment of 0 should always be 10

PCprey     = EVprey./sum(EVprey);                       % compute EV of the offers given current model
PCpred     = EVpred./sum(EVpred);     


figure
set(gcf,'Color',[1,1,1])
subplot(2,2,1)
hold on
plot(offers,XPred,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,XPrey,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('opponent probability')
legend('prey','predator')
title('Level 2')

subplot(2,2,2)
hold on
plot(offers,YPred,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,YPrey,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('p(X>(opponent))')
legend('prey','predator')

subplot(2,2,3)
hold on
plot(offers,EVprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,EVpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('EV')
legend('prey','predator')

subplot(2,2,4)
hold on
plot(offers,PCprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,PCpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('likelihood')
legend('prey','predator')


%% level-3
% at level 3, players expect opponent to be level-2 players

XPrey = PCprey;
XPred = PCpred;

YPrey = cumsum(XPrey);
YPred = cumsum(XPred);

EVprey     = (endow - offers).* (YPred);                       % compute EV of the offers given current model
% EVpred     = (endow - offers) + (endow - offers).* (Y);   % compute EV of the offers given current model
EVpred = NaN(1, length(offers));
for tc = 1:length(EVpred)
    EVpred(tc) = ((endow(tc) - offers(tc))) + sum((endow(1:tc) - offers(1:tc)) .* (XPred(1:tc)));
end
EVpred(1) = 10; % EV for an investment of 0 should always be 10

PCprey     = EVprey./sum(EVprey);                       % compute EV of the offers given current model
PCpred     = EVpred./sum(EVpred);     


figure
set(gcf,'Color',[1,1,1])
subplot(2,2,1)
hold on
plot(offers,XPred,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,XPrey,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('opponent probability')
legend('prey','predator')
title('Level 3')

subplot(2,2,2)
hold on
plot(offers,YPred,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,YPrey,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('p(X>(opponent))')
legend('prey','predator')

subplot(2,2,3)
hold on
plot(offers,EVprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,EVpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('EV')
legend('prey','predator')

subplot(2,2,4)
hold on
plot(offers,PCprey,'-o',...
    'MarkerFaceColor',[1,1,1])
plot(offers,PCpred,'-o',...
    'MarkerFaceColor',[1,1,1])
xlabel('offers')
ylabel('likelihood')
legend('prey','predator')
