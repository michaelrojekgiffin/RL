clear
close all force
clc


%% load priors
load('Learning_All_2017_12_08')

%% find paths
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');
fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));
nsub = length(fl_list);



%% some params
nfpm = [4,5,4,5];
Nmodel = 2;
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

%% model comparison
X1 = squeeze(LEARN_LPP(:,:,1));
for k_mod= 1:4
    bic(:,k_mod)=-2*-X1(:,k_mod) + nfpm(k_mod)*log(3*25); % l2 is already positive
end
[postBMC1,outBMC1]=VBA_groupBMC(-bic'./2);

X2 = squeeze(LEARN_LPP(:,:,2));
for k_mod= 1:4
    bic(:,k_mod)=-2*-X2(:,k_mod) + nfpm(k_mod)*log(3*25); % l2 is already positive
end
[postBMC2,outBMC2]=VBA_groupBMC(-bic'./2);

count1 = 0;
for k_sub = 1:nsub;
    
    %% check sub
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        %% Get UG data
        flnm = fullfile(data_dir,strcat('Sub_',sub_nm,'_UG.mat'));
        load(flnm)
        impexp(k_sub) = subdata(1,6);
        
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
            
            [PE,eO,a_t] = learning_models_outputs(LEARN_parametersLPP(k_sub,:,Nmodel,kSoc),O,D,Nmodel);
            
            eO_mat(k_sub,:,:,kSoc) = eO;
            logA_mat(k_sub,:,kSoc) = a_t(end,:);
            logBmat(k_sub,:,kSoc) = repmat(LEARN_parametersLPP(k_sub,3,Nmodel,kSoc),1,3);
            
            for kcond = 1:3
                PA_sub(k_sub,:,kcond,kSoc)     = logitp([squeeze(logA_mat(k_sub,kcond,kSoc))',squeeze(logBmat(k_sub,kcond,kSoc))'],offers);            % compute proba of accepting the offers given current model
                EV_sub(k_sub,:,kcond,kSoc)    = (endow - offers).* squeeze(PA_sub(k_sub,:,kcond,kSoc)) ;
            end
        end
    end
end




colmat = [1,0,0;.5,0,.5;0,0,1];
symb_IE = {'o','sq'};

%% Plot learning curves

figure;
set(gcf,'Color',[1,1,1])
for kSoc = 1:2
    
    for kImpexp = 1:2
        
        subplot(2,2,(kSoc-1)*2 + kImpexp)
        
        Gsub = impexp == (kImpexp -1);
        nGsub = sum(double(Gsub));
        
        easy_sub{kSoc,kImpexp} = LEARN_parametersLPP(Gsub,:,Nmodel,kSoc);
        
        hold on
        bar(squeeze(mean(LEARN_parametersLPP(Gsub,:,Nmodel,kSoc),1)),...
            'FaceColor',.7*[1,1,1])
        errorbar(squeeze(mean(LEARN_parametersLPP(Gsub,:,Nmodel,kSoc),1)),squeeze(std(LEARN_parametersLPP(Gsub,:,Nmodel,kSoc),0,1))./sqrt(nGsub),'k',...
            'LineStyle','none')
        
        set(gca,'XLim',[0 6],...
            'XTick',1:5,...
            'XTickLabel',{'\beta','a_0','b_0','\alpha_1','\alpha_2'},...
            'YLim',[-5 5])
        xlabel('param')
        ylabel('Value')
        title(strcat(['Social = ',num2str(kSoc-1) '; ImpExp == ',num2str(kImpexp-1)]));
    end
end

% pNSoc = squeeze(LEARN_parametersLPP(Gsub,:,Nmodel,1));
% pSoc = squeeze(LEARN_parametersLPP(Gsub,:,Nmodel,2));

pNSoc = [easy_sub{1,1};easy_sub{1,2}];
pSoc  = [easy_sub{2,1};easy_sub{2,2}];
[H P CI STATS] = ttest(pNSoc,pSoc)

pImp = [easy_sub{1,1}+easy_sub{2,1}]./2;
pExp  = [easy_sub{1,2}+easy_sub{2,2}]./2;
[H P CI STATS] = ttest2(pImp,pExp)

%% Plot Average Params


figure;
set(gcf,'Color',[1,1,1])


for kSoc = 1:2
    
    for kImpexp = 1:2
        
        subplot(2,2,(kSoc-1)*2 + kImpexp)
        Gsub = impexp == (kImpexp -1);
        nGsub = sum(double(Gsub));
        
        hold on
        for k_cond = 1:3
            
            Ymod = squeeze(eO_mat(Gsub,:,k_cond,kSoc));
            fill([(1:ntr),fliplr(1:ntr)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nGsub)))],...
                colmat(k_cond,:),'LineStyle','none')
            alpha(.15)
            %        plot(mean(Ymod),'-',...
            %    'Color',colmat(k_cond,:))
            
            errorbar(squeeze(mean(offer_mat(Gsub,:,k_cond,kSoc),1)),squeeze(std(offer_mat(Gsub,:,k_cond,kSoc),0,1))./sqrt(nGsub),symb_IE{kImpexp},...
                'Color',colmat(k_cond,:),...
                'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
                'MarkerEdgeColor',colmat(k_cond,:))
        end
        set(gca,'XLim',[0 ntr+1],...
            'YLim',[4 10])
        xlabel('Trial')
        ylabel('Offer')
        %  legend('Cond1','Cond2','Cond3')
        title(strcat(['Social = ',num2str(kSoc-1) '; ImpExp == ',num2str(kImpexp-1)]));
    end
end

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
        
        % plot final model acceptance estimations
        Ymod = squeeze(PA_sub(Gsub,:,k_cond,kSoc));
        fill([(offers),fliplr(offers)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nsub)))],...
            colmat(k_cond,:),'LineStyle','none')
        alpha(.15)
        plot(offers,mean(Ymod),'-',...
            'Color',colmat(k_cond,:))
        
        % plot observed acceptance frequency (~estimation if subjects are optimal)
        plot(nOFF,off_rate,'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
        
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
        
        % plot model
        Ymod = squeeze(EV_sub(Gsub,:,k_cond,kSoc));
        plot(offers,mean(Ymod),'-',...
            'Color',colmat(k_cond,:))
        fill([(offers),fliplr(offers)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nsub)))],...
            colmat(k_cond,:),'LineStyle','none')
        alpha(.15)
        
        plot(nOFF,off_rate,'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
        
        
        
        set(gca,'YLim',[0 1])
    end
    
    set(gca,'XLim',[0 20],...
        'YLim',[0 20])
    xlabel('Offer')
    ylabel('Expected value')
    title(strcat(['Social = ',num2str(kSoc-1)]))
    
end



