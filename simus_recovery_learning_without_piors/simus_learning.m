clear
close all
clc

% randomize generator seed
%--------------------------
rng('shuffle')

% parameters of the task
%--------------------------
n_sims      = 100;                           % nsubs to simulates
n_trial     = 24;                           % ntrial per cond per session
n_sess      = 100;                            % nsession
offers      = 0:1:20;
suboffers   = 0:.001:20;
endow       = 20*ones(1,numel(offers));% parameters of the simulation
subendow    = 20*ones(1,numel(suboffers));% parameters of the simulation
nmodel      = 2;

% set up conditions and mutliple sessions
%------------------------------------------
cond2learn  = -[9,6,3];
% cond2learn  = -[10,8,6,4,2];
nc          = numel(cond2learn);
Ra          = repmat(cond2learn,1,n_sess);            % predator true accepance thereshold (logit intercept)
Rb          = repmat(1*ones(1,nc),1,n_sess);              % predator true accpetance moise (logit slope)
n_cond      = size(Ra,2);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Generate params
%-------------------
Pa_rnd          = round(100*linspace(-12,0,n_sims+1))./100;  %  Proposer initial prior on threshold (to be fitted or estimated from no-feedback games)
Pb_rnd          = 1*ones(n_sims,1);       %  Proposer  estimated accpetance noise (to be fitted or estimated from no-feedback games)
Px_rnd          = 5*ones(n_sims,1);         %  Proposer  rating temperature
Plr1_rnd        = .25*ones(n_sims,1);           %  Proposer  learning rate
Plr2_rnd        = .1*ones(n_sims,1);           %  Proposer  learning rate


% PreAllocate
%---------------
As_mat = zeros(n_sims,n_trial,nc);
A_mat = zeros(n_sims,n_trial,nc);
PE_mat = zeros(n_sims,n_trial,nc);
O_mat = zeros(n_sims,n_trial,nc);
% Sim loop
%----------
for k_sim = 1:n_sims
    
    % get params
    a0  = Pa_rnd(k_sim);
    b0  = Pb_rnd(k_sim);
    bX  = Px_rnd(k_sim);
    lr1 = Plr1_rnd(k_sim);
    lr2 = Plr2_rnd(k_sim);
    
    [O,D,PE,at] = learning_models_timeseries([bX,a0,b0,lr1,lr2],[Ra;Rb],n_trial,nmodel);
    A_mat(k_sim,:,:) = squeeze(mean(reshape(at(1:n_trial,:),n_trial,nc,n_sess),3));
    PE_mat(k_sim,:,:) = squeeze(mean(reshape(PE(1:n_trial,:),n_trial,nc,n_sess),3));
    O_mat(k_sim,:,:) = squeeze(mean(reshape(O(1:n_trial,:),n_trial,nc,n_sess),3));
    As_mat(k_sim,:,:) = squeeze(std(reshape(at(1:n_trial,:),n_trial,nc,n_sess),0,3))./sqrt(n_sess);
end

save('Simus_2018_03_12_3')

%%
yleg = [-12:3:0];
for kt = 1:numel(yleg)
    [~,yline(kt)] = min(abs(Pa_rnd-yleg(kt)));
end

h1 = figure('Units', 'pixels', ...
    'Position', [400 300 1200 250]);
set(gcf,'Color',[1,1,1])
%[map] = cbrewer('div', 'Spectral',11);
[map] = cbrewer('div', 'BrBG',11);
colormap(map)

for k = 1:nc
    subplot(1,nc,k)
    mtp = squeeze(A_mat(:,:,k));
    
    [X,Y] = meshgrid(1:1:n_trial,1:1:n_sims);
    Z = (mtp - cond2learn(k).*ones(n_sims,n_trial));
    contourf(X,Y,Z,11,'LineStyle','none')
    
    hX = xlabel('Trials');
    hY = ylabel('Initial a');
    hT = title(strcat('True a = ',num2str(cond2learn(k))));
    set(gca,'YTick',yline,...
        'YTickLabel',yleg,...
        'FontName','Arial',...
        'FontSize',9)
    set([hX,hY],'FontName','Arial',...
        'FontSize',9)
    
    colorbar
    caxis([-10,10])
end
% print('Test','-dsvg')
  print(gcf,'-depsc','-painters','out.eps');
       epsclean('out.eps'); % cleans and overwrites the input file

h2 = figure('Units', 'pixels', ...
    'Position', [400 300 1200 250]);
set(gcf,'Color',[1,1,1])
[map] = cbrewer('seq', 'BuGn',15);
colormap(map)
for k = 1:nc
    subplot(1,nc,k)
    mtp = squeeze(O_mat(:,:,k));
    [~,m] = max((subendow - suboffers).*logitp([cond2learn(k),1],suboffers));             % compute EV of the offers given current model
    best_o = suboffers(m);
    [X,Y] = meshgrid(1:1:n_trial,1:1:n_sims);
    Z = mtp;
    contourf(X,Y,Z,11,'LineStyle','none')
    
    hX = xlabel('Trials');
    hY = ylabel('Initial Offer');
    hT = title(strcat('bst Offer = ',num2str(best_o)));
    set(gca,'YTick',yline,...
        'YTickLabel',yleg,...
        'FontName','Arial',...
        'FontSize',9)
    set([hX,hY],'FontName','Arial',...
        'FontSize',9)
    
    colorbar
    caxis([0,15])
end

% h1 = figure('Units', 'pixels', ...
%     'Position', [400 300 1200 250]);
% set(gcf,'Color',[1,1,1])

% for k = 1:nc
%     subplot(2,nc,nc+k)
%     mtp = squeeze(O_mat(:,:,k));
%     [~,m] = max((subendow - suboffers).*logitp([cond2learn(k),1],suboffers));             % compute EV of the offers given current model
%     best_o = suboffers(m);
%     [X,Y] = meshgrid(1:1:n_trial,1:1:n_sims);
%     Z = (mtp - best_o.*ones(n_sims,n_trial));
%     contourf(X,Y,Z,11,'LineStyle','none')
%
%     hX = xlabel('Trials');
%     hY = ylabel('Initial a');
%     hT = title(strcat('True a = ',num2str(cond2learn(k))));
%     set(gca,'YTick',yline,...
%         'YTickLabel',yleg,...
%         'FontName','Arial',...
%         'FontSize',9)
%     set([hX,hY],'FontName','Arial',...
%         'FontSize',9)
%
%     colorbar
%     caxis([-5,5])
% end

%%
h3 = figure('Units', 'pixels', ...
    'Position', [400 300 1200 250]);
set(gcf,'Color',[1,1,1])

for k = 1:nc
    Cmap(k,:) = [(k-1)/(nc-1),0,(nc-k)./(nc-1)];
end
colormap(Cmap)    
    
condstp = [1 25 50];
for kc = 1:numel(condstp)
    
    subplot(1,3,kc)
    hold on
    ctp = condstp(kc);
    for k = 1:nc
        errorbar(squeeze(A_mat(ctp,:,k)),squeeze(As_mat(ctp,:,k)),'-o',...
            'MarkerFaceColor',Cmap(k,:),...
            'MarkerEdgeColor',Cmap(k,:),...
            'Color',Cmap(k,:));
        lg{k} = strcat(num2str(cond2learn(k))); 
    end
    
    hX = xlabel('Trials');
    hY = ylabel('Initial a');
    set(gca,'YLim',[-12 0],...
        'XLim',[0 n_trial+1],...
        'FontName','Arial',...
        'FontSize',9)
    
    colorbar
    caxis([1,nc])
    colorbar('Ticks',1:nc,...
        'TickLabels',lg)
end


