clear
close all force
clc

%% find paths
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');
fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));
nsub = length(fl_list);

%% load priors
load('Priors_All_2017_12_07')
load('Learning_All_2017_12_07')

%% some params
Nmodel = 1;
ntr = 25;
ncond = 3;

%% pre allocate
good_sub = NaN(nsub,1);

offer_mat = NaN(nsub,ntr,ncond,2);
accept_mat = NaN(nsub,ntr,ncond,2);
reward_mat = NaN(nsub,ntr,ncond,2);
cond_mat = NaN(nsub,ntr,ncond,2);
sub_mat = NaN(nsub,ntr,ncond,2);
soc_mat = NaN(nsub,ntr,ncond,2);
trial_mat = NaN(nsub,ntr,ncond,2);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:20;
endow  = 20*ones(1,numel(offers));% parameters of the simulation



for k_sub = 1:nGsub;
    
    %% check sub
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        good_sub(k_sub) = 1;
        %% Get UG data
        flnm = fullfile(data_dir,strcat('Sub_',sub_nm,'_UG.mat'));
        load(flnm)
        
        for kSoc = 1:2
            SOC = subdata(subdata(:,4) == kSoc-1,:);
            cond1 = SOC(SOC(:,5)==1,:);
            cond2 = SOC(SOC(:,5)==2,:);
            cond3 = SOC(SOC(:,5)==3,:);
            offer_mat(k_sub,:,:,kSoc) = [cond1(1:ntr,1),cond2(1:ntr,1),cond3(1:ntr,1)];
            accept_mat(k_sub,:,:,kSoc) = [cond1(1:ntr,2),cond2(1:ntr,2),cond3(1:ntr,2)];
            reward_mat(k_sub,:,:,kSoc) = [cond1(1:ntr,3),cond2(1:ntr,3),cond3(1:ntr,3)];
            cond_mat(k_sub,:,:,kSoc) = [ones(ntr,1),2*ones(ntr,1),3*ones(ntr,1)];
            sub_mat(k_sub,:,:,kSoc) = k_sub*ones(ntr,3);
            soc_mat(k_sub,:,:,kSoc) = (kSoc)*ones(ntr,3);
            trial_mat(k_sub,:,:,kSoc) = repmat((1:ntr)',1,3);
            
            
            O = squeeze(offer_mat(k_sub,:,:,kSoc));
            D = squeeze(accept_mat(k_sub,:,:,kSoc));
            a0 = squeeze(parametersLPP(k_sub,1,kSoc));
            b0 = squeeze(parametersLPP(k_sub,2,kSoc));
            
            [PE,eO,a_t] = learning_models_outputs(LEARN_parametersLPP(k_sub,:,Nmodel,kSoc),O,D,a0,b0,Nmodel);
            
            eO_mat(k_sub,:,:,kSoc) = eO;
            logA_mat(k_sub,:,kSoc) = a_t(end,:);
            logBmat(k_sub,:,kSoc) = repmat(b0,1,3);
            
            for kcond = 1:3
                PA_sub(k_sub,:,kcond,kSoc)     = logitp([squeeze(logA_mat(k_sub,kcond,kSoc))',squeeze(logBmat(k_sub,kcond,kSoc))'],offers);            % compute proba of accepting the offers given current model
                EV_sub(k_sub,:,kcond,kSoc)    = (endow - offers).* squeeze(PA_sub(k_sub,:,kcond,kSoc)) ;
            end
        end
    end
end





Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

colmat = [1,0,0;.5,0,.5;0,0,1];

%% Plot learning curves

figure;
set(gcf,'Color',[1,1,1])


for kSoc = 1:2
    
    subplot(1,2,kSoc)
    hold on
    for k_cond = 1:3
        errorbar(squeeze(mean(offer_mat(Gsub,:,k_cond,kSoc),1)),squeeze(std(offer_mat(Gsub,:,k_cond,kSoc),0,1))./sqrt(nGsub),'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
        
        plot(squeeze(mean(eO_mat(Gsub,:,k_cond,kSoc),1)),'-',...
            'Color',colmat(k_cond,:))
    end
    set(gca,'XLim',[0 ntr+1],...
        'YLim',[4 10])
    xlabel('Trial')
    ylabel('Offer')
    legend('Cond1','Cond2','Cond3')
    title(strcat(['Social = ',num2str(kSoc-1)]));
end


%% Plot Average offers

X_non_SOC = squeeze(mean(offer_mat(Gsub,:,:,1),2));
X_SOC = squeeze(mean(offer_mat(Gsub,:,:,2),2));

mtp_SOC = squeeze(mean(X_SOC,1));
stp_SOC = squeeze(std(X_SOC,0,1))./sqrt(nGsub);

mtp_non_SOC = squeeze(mean(X_non_SOC,1));
stp_non_SOC = squeeze(std(X_non_SOC,0,1))./sqrt(nGsub);

figure;
set(gcf,'Color',[1,1,1])
subplot(1,2,1)
hold on

for k_cond = 1:3
    bar(k_cond,mtp_non_SOC(:,k_cond),...
        'FaceColor',[1,1,1],...
        'EdgeColor',colmat(k_cond,:))
    errorbar(k_cond,mtp_non_SOC(:,k_cond),stp_non_SOC(:,k_cond),'k','LineStyle','none')
end

set(gca,'XTick',1:3,...
    'YLim',[4 10])
xlabel('COND')
ylabel('Average Offer')
title('Social = 0')

subplot(1,2,2)
hold on

for k_cond = 1:3
    bar(k_cond,mtp_SOC(:,k_cond),...
        'FaceColor',colmat(k_cond,:),...
        'EdgeColor',[0,0,0])
    errorbar(k_cond,mtp_SOC(:,k_cond),stp_SOC(:,k_cond),'k','LineStyle','none')
end
set(gca,'XTick',1:3,...
    'YLim',[4 10])
xlabel('COND')
ylabel('Average Offer')
title('Social = 1')


%% Plot observed acceptance rate
figure
set(gcf,'Color',[1,1,1])

for kSoc = 1:2
    
    subplot(1,2,kSoc)
    hold on
    
    for k_cond = 1:3
        
        offer = squeeze(offer_mat(Gsub,:,k_cond,kSoc));
        offer = offer(:);
        accept = squeeze(accept_mat(Gsub,:,k_cond,kSoc));
        accept = accept(:);
        
        nOFF = unique(offer);
        off_rate =[];
        for koff = 1:length(nOFF);
            off_rate(koff) = mean(accept(offer == koff));
        end
        
        plot(nOFF,off_rate,'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
        
        plot(offers,squeeze(mean(PA_sub(Gsub,:,k_cond,kSoc))),'-',...
            'Color',colmat(k_cond,:))
        
        set(gca,'YLim',[0 1])
    end
    
    set(gca,'XLim',[0 20],...
        'YLim',[0 1])
    xlabel('Offer')
    ylabel('Acceptance rate')
    title(strcat(['Social = ',num2str(kSoc-1)]))
    
end

%% Plot expected reward
figure
set(gcf,'Color',[1,1,1])

for kSoc = 1:2
    
    subplot(1,2,kSoc)
    hold on
    
    for k_cond = 1:3
        
        offer = squeeze(offer_mat(Gsub,:,k_cond,kSoc));
        offer = offer(:);
        reward = squeeze(reward_mat(Gsub,:,k_cond,kSoc));
        reward = reward(:);
        
        nOFF = unique(offer);
        off_rate =[];
        for koff = 1:length(nOFF);
            off_rate(koff) = mean(reward(offer == koff));
        end
        
        plot(nOFF,off_rate,'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
        
        plot(offers,squeeze(mean(EV_sub(Gsub,:,k_cond,kSoc))),'-',...
            'Color',colmat(k_cond,:))
        
        set(gca,'YLim',[0 1])
    end
    
    set(gca,'XLim',[0 20],...
        'YLim',[0 20])
    xlabel('Offer')
    ylabel('Expected value')
    title(strcat(['Social = ',num2str(kSoc-1)]))
    
end



