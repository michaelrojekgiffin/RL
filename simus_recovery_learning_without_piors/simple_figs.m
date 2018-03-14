clear all
close all
clc



offers      = 0:1:20;
suboffers   = 0:.1:20;
endow       = 20*ones(1,numel(offers));% parameters of the simulation
subendow    = 20*ones(1,numel(suboffers));% parameters of the simulation

P = [-6 1];

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

logitd = @(b,x) exp(-(x-b(1))./b(2))./(b(2).*(1+exp(-(x-b(1))./b(2)).^2));

% Generate params

V = (subendow - suboffers);
PA = logitp(P,suboffers);
EV = (subendow - suboffers).*logitp(P,suboffers);


h1 = figure('Units', 'pixels', ...
    'Position', [400 300 800 200]);
set(gcf,'Color',[1,1,1])

subplot(1,3,1)
hold on
plot(suboffers,V,'-')
plot(suboffers,zeros(1,numel(suboffers)),'-')
hX = xlabel('Offer');
hY = ylabel('Value');
set(gca,'XLim',[0 20],...
    'YLim',[0 20],...
    'FontName','Arial',...
    'FontSize',9)

subplot(1,3,2)
plot(suboffers,100*PA,'-')
hX = xlabel('Offer');
hY = ylabel('P(accepted) ');
set(gca,'XLim',[0 20],...
    'YLim',[0 100],...
    'FontName','Arial',...
    'FontSize',9)

subplot(1,3,3)
plot(suboffers,EV,'-')
hX = xlabel('Offer');
hY = ylabel('Expected Value');
set(gca,'XLim',[0 20],...
    'YLim',[0 12],...
    'FontName','Arial',...
    'FontSize',9)


% toPPT(h1,'format','vec')



%%
h2a = figure('Units', 'pixels', ...
    'Position', [400 300 250 200]);
set(gcf,'Color',[1,1,1])

hold on
PA1 = logitp([-9 1],suboffers);
PA2 = logitp([-6 1],suboffers);
PA3 = logitp([-3 1],suboffers);
plot(suboffers,100*PA1,'-')
plot(suboffers,100*PA2,'-')
plot(suboffers,100*PA3,'-')
hX = xlabel('Offer');
hY = ylabel('P(accepted) ');
set(gca,'XLim',[0 20],...
    'YLim',[0 100],...
    'FontName','Arial',...
    'FontSize',9)
legend('-9','-6','-3')


%%
h2b = figure('Units', 'pixels', ...
    'Position', [400 300 250 200]);
set(gcf,'Color',[1,1,1])

hold on
PA1 = logitp([-9 1],suboffers);
PA2 = logitp([-6 1],suboffers);
PA3 = logitp([-3 1],suboffers);
plot(suboffers,100*PA1,'-','Color',[0,0,1])
plot(suboffers,100*PA2,'-','Color',[.5,0,.5])
plot(suboffers,100*PA3,'-','Color',[1,0,0])
hX = xlabel('Offer');
hY = ylabel('P(accepted) ');
set(gca,'XLim',[0 20],...
    'YLim',[0 100],...
    'FontName','Arial',...
    'FontSize',9)
legend('-9','-6','-3')




% toPPT(h2b,'format','vec')


h3 = figure('Units', 'pixels', ...
    'Position', [400 300 250 200]);
set(gcf,'Color',[1,1,1])

hold on
PA1 = logitd([9./1 1],suboffers);
PA2 = logitd([6 1],suboffers);
PA3 = logitd([3 1],suboffers);
plot(suboffers,100*PA1,'-')
plot(suboffers,100*PA2,'-')
plot(suboffers,100*PA3,'-')
hX = xlabel('Offer');
hY = ylabel('P(accepted) ');
set(gca,'XLim',[0 20],...
    'YLim',[0 80],...
    'FontName','Arial',...
    'FontSize',9)
legend('-9','-6','-3')

% toPPT(h3,'format','vec')
%%
h4 = figure('Units', 'pixels', ...
    'Position', [400 300 500 200]);
set(gcf,'Color',[1,1,1])

pEV = (endow - offers).*logitp(P,offers);
subplot(1,2,1)
plot(offers,pEV,'-')
hX = xlabel('Offer');
hY = ylabel('Expected Value');
set(gca,'XLim',[0 20],...
    'YLim',[0 12],...
    'FontName','Arial',...
    'FontSize',9)

subplot(1,2,2)
hold on
PC1 = exp(.25.*pEV) ./ sum(exp(.25.*pEV));   %
PC2 = exp(1.*pEV) ./ sum(exp(1.*pEV));   %
PC3 = exp(5.*pEV) ./ sum(exp(5.*pEV));   %
plot(offers,100*PC1,'-')
plot(offers,100*PC2,'-')
plot(offers,100*PC3,'-')
hX = xlabel('Offer');
hY = ylabel('probability choice');
set(gca,'XLim',[0 20],...
    'YLim',[0 80],...
    'FontName','Arial',...
    'FontSize',9)
legend('0.25','1','5')



h5 = figure('Units', 'pixels', ...
    'Position', [400 300 500 200]);
set(gcf,'Color',[1,1,1])

subplot(1,2,1)
plot(suboffers,EV,'-')
hX = xlabel('Offer');
hY = ylabel('Expected Value');
set(gca,'XLim',[0 20],...
    'YLim',[0 12],...
    'FontName','Arial',...
    'FontSize',9)

subplot(1,2,2)
hold on
PC1 = exp(.25.*EV) ./ sum(exp(.25.*EV));   %
PC2 = exp(1.*EV) ./ sum(exp(1.*EV));   %
PC3 = exp(5.*EV) ./ sum(exp(5.*EV));   %
plot(suboffers,100*PC1,'-')
plot(suboffers,100*PC2,'-')
plot(suboffers,100*PC3,'-')
hX = xlabel('Offer');
hY = ylabel('probability choice');
set(gca,'XLim',[0 20],...
    'YLim',[0 10],...
    'FontName','Arial',...
    'FontSize',9)
legend('0.25','1','5')



 % toPPT(h5,'format','vec')