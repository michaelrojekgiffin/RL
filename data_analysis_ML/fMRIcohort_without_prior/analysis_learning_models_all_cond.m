clear
close all force
clc


%% load priors
load('Learning_All_2018_02_15')

%% find paths
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'fMRI_stimulus','UG','data');
subjects = [3:14];
nsub = numel(subjects);



%% some params
nfpm = [4,5,4,5];
Nmodel = 2;
ncond = 6;

%% pre allocate
good_sub = NaN(nsub,1);
ntr = 24;
offer_mat_sess = NaN(nsub,ntr,ncond,2);
accept_mat_sess = NaN(nsub,ntr,ncond,2);
reward_mat_sess = NaN(nsub,ntr,ncond,2);
cond_mat_sess = NaN(nsub,ntr,ncond,2);
sub_mat_sess = NaN(nsub,ntr,ncond,2);
soc_mat_sess = NaN(nsub,ntr,ncond,2);
trial_mat_sess = NaN(nsub,ntr,ncond,2);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:20;
endow  = 20*ones(1,numel(offers));% parameters of the simulation

%% model comparison
X1 = squeeze(LEARN_LPP(:,:,1));
for k_mod= 1:4
    bic(:,k_mod)=-2*-X1(:,k_mod) + nfpm(k_mod)*log(6*24); % l2 is already positive
end
[postBMC1,outBMC1]=VBA_groupBMC(-bic'./2);

X2 = squeeze(LEARN_LPP(:,:,2));
for k_mod= 1:4
    bic(:,k_mod)=-2*-X2(:,k_mod) + nfpm(k_mod)*log(6*24); % l2 is already positive
end
[postBMC2,outBMC2]=VBA_groupBMC(-bic'./2);

count1 = 0;
%% Loop Subjects
for k_sub = 1:nsub;
    
    
    full_data = NaN(2*ntr,7);
    sub_num = sprintf('%03.0f',(subjects(k_sub)));
    
    for kSoc = 1:2
        
        for k_sess = 1:2
            k_in = (k_sess-1)*ntr+1;
            k_out = k_sess*ntr;
            %% Get UG data
            switch kSoc
                case 2
                    fldir = dir(strcat(data_dir,filesep,'sub',sub_num,'_social_*'));
                case 1
                    fldir = dir(strcat(data_dir,filesep,'sub',sub_num,'_nonsoc_*'));
            end
            
            flnm = fullfile(data_dir,fldir(k_sess).name);
            if strcmp(sub_num,'013') && k_sess == 2 && kSoc ==2
                flnm = fullfile(data_dir,fldir(3).name);
            end
            load(flnm)
            
            cond1 = sub_data(sub_data(:,5)==0,:);
            cond2 = sub_data(sub_data(:,5)==1,:);
            cond3 = sub_data(sub_data(:,5)==2,:);
            
            switch k_sess
                case 1
                    offer_mat_sess(k_sub,:,1:3,kSoc) = [cond1(:,1),cond2(:,1),cond3(:,1)];
                    accept_mat_sess(k_sub,:,1:3,kSoc) = [cond1(:,2),cond2(:,2),cond3(:,2)];
                    reward_mat_sess(k_sub,:,1:3,kSoc) = [cond1(:,3),cond2(:,3),cond3(:,3)];
                    cond_mat_sess(k_sub,:,1:3,kSoc) = [ones(24,1),2*ones(24,1),3*ones(24,1)];
                    sub_mat_sess(k_sub,:,1:3,kSoc) = k_sub*ones(24,3);
                    soc_mat_sess(k_sub,:,1:3,kSoc) = (kSoc-1)*ones(24,3);
                    trial_mat_sess(k_sub,:,1:3,kSoc) = repmat((1:24)',1,3);
                case 2
                    offer_mat_sess(k_sub,:,4:6,kSoc) = [cond1(:,1),cond2(:,1),cond3(:,1)];
                    accept_mat_sess(k_sub,:,4:6,kSoc) = [cond1(:,2),cond2(:,2),cond3(:,2)];
                    reward_mat_sess(k_sub,:,4:6,kSoc) = [cond1(:,3),cond2(:,3),cond3(:,3)];
                    cond_mat_sess(k_sub,:,4:6,kSoc) = 3+[ones(24,1),2*ones(24,1),3*ones(24,1)];
                    sub_mat_sess(k_sub,:,4:6,kSoc) = k_sub*ones(24,3);
                    soc_mat_sess(k_sub,:,4:6,kSoc) = (kSoc-1)*ones(24,3);
                    trial_mat_sess(k_sub,:,4:6,kSoc) = repmat((1:24)',1,3);
                    
            end
            
        end
        
        O = squeeze(offer_mat(k_sub,:,:,kSoc));
        D = squeeze(accept_mat(k_sub,:,:,kSoc));
        
        
        
        [PE,eO,a_t] = learning_models_outputs(LEARN_parametersLPP(k_sub,:,Nmodel,kSoc),O,D,Nmodel);
        
        eO_mat_sess(k_sub,:,:,kSoc) = eO;
        logA_mat_sess(k_sub,:,kSoc) = a_t(end,:);
        logB_mat_sess(k_sub,:,kSoc) = repmat(LEARN_parametersLPP(k_sub,3,Nmodel,kSoc),1,6);
        
        for kcond = 1:6
            PA_sub_sess(k_sub,:,kcond,kSoc)     = logitp([squeeze(logA_mat_sess(k_sub,kcond,kSoc))',squeeze(logB_mat_sess(k_sub,kcond,kSoc))'],offers);            % compute proba of accepting the offers given current model
            EV_sub_sess(k_sub,:,kcond,kSoc)    = (endow - offers).* squeeze(PA_sub_sess(k_sub,:,kcond,kSoc)) ;
        end
        
    end
end

eO_mat = (eO_mat_sess(:,:,1:3,:) + eO_mat_sess(:,:,4:6,:))./2;
logA_mat = (logA_mat_sess(:,1:3,:) + logA_mat_sess(:,4:6,:))./2;
logB_mat = (logA_mat_sess(:,1:3,:) + logA_mat_sess(:,4:6,:))./2;

PA_sub =(PA_sub_sess(:,:,1:3,:)+PA_sub_sess(:,:,4:6,:))./2;
EV_sub =(EV_sub_sess(:,:,1:3,:)+EV_sub_sess(:,:,4:6,:))./2;



offer_mat = (offer_mat_sess(:,:,1:3,:) + offer_mat_sess(:,:,4:6,:))./2;
accept_mat = (accept_mat_sess(:,:,1:3,:) + accept_mat_sess(:,:,4:6,:))./2;
reward_mat = (reward_mat_sess(:,:,1:3,:) + reward_mat_sess(:,:,4:6,:))./2;
cond_mat = (cond_mat_sess(:,:,1:3,:) + cond_mat_sess(:,:,4:6,:))./2;
sub_mat = (sub_mat_sess(:,:,1:3,:) + sub_mat_sess(:,:,4:6,:))./2;
soc_mat = (soc_mat_sess(:,:,1:3,:) + soc_mat_sess(:,:,4:6,:))./2;
trial_mat = (trial_mat_sess(:,:,1:3,:) + trial_mat_sess(:,:,4:6,:))./2;


%%

colmat = [1,0,0;.5,0,.5;0,0,1];

%% Plot params

h1 = figure('Units', 'pixels', ...
    'Position', [400 200 800 250]);
set(gcf,'Color',[1,1,1])

yL = [0 15;...
    -12 0;...
    0 3;...
    0 1;...
    0 0.2];

Pnames = {'\beta','a_0','b_0','\alpha_1','\alpha_2'};

for k_param = 1:5
    
    
    subplot(1,5,k_param)
    hold on
    for k_sub = 1:nsub
        plot([.8 2.2],[LEARN_parametersLPP(k_sub,k_param,Nmodel,2) LEARN_parametersLPP(k_sub,k_param,Nmodel,1)],...
            '-k','LineWidth',.5,'Color',.7*[1,1,1])
    end
    plot(.8*ones(nsub,1),squeeze(LEARN_parametersLPP(:,k_param,Nmodel,2)),'o',...
        'LineStyle','none',...
        'MarkerEdgeColor',.5.*[1,1,1],...
        'MarkerFaceColor',.6.*[1,1,1])
    plot(2.2*ones(nsub,1),squeeze(LEARN_parametersLPP(:,k_param,Nmodel,1)),'o',...
        'LineStyle','none',...
        'MarkerEdgeColor',.7.*[1,1,1],...
        'MarkerFaceColor',[1,1,1])
    
    errorbar([1.25 1.75],...
        [squeeze(mean(LEARN_parametersLPP(:,k_param,Nmodel,2),1)) squeeze(mean(LEARN_parametersLPP(:,k_param,Nmodel,1),1))],...
        [squeeze(std(LEARN_parametersLPP(:,k_param,Nmodel,2),0,1))./sqrt(nsub) squeeze(std(LEARN_parametersLPP(:,k_param,Nmodel,1),0,1))./sqrt(nsub)],...
        '-k','LineWidth',1,'Color',.2*[1,1,1])
    plot(1.25,squeeze(mean(LEARN_parametersLPP(:,k_param,Nmodel,2),1)),'o',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',.0.*[1,1,1])
    plot(1.75,squeeze(mean(LEARN_parametersLPP(:,k_param,Nmodel,1),1)),'o',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',[1,1,1])
  
    
    hT = title(Pnames{k_param});
    
    set(gca,'XLim',[0.25 2.75],...
        'YLim',yL(k_param,:),...
        'XTick', [1 2],...
        'XTickLabel',{'S','nS'})
end



%% Plot learning curves


figure;
set(gcf,'Color',[1,1,1])


for kSoc = 1:2
    
    subplot(1,2,kSoc)
    Gsub = ~isnan(subjects);
    nGsub = sum(double(Gsub));
    
    hold on
    for k_cond = 1:3
        
        Ymod = squeeze(eO_mat(Gsub,:,k_cond,kSoc));
        fill([(1:ntr),fliplr(1:ntr)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nGsub)))],...
            colmat(k_cond,:),'LineStyle','none')
        alpha(.15)
        %        plot(mean(Ymod),'-',...
        %    'Color',colmat(k_cond,:))
        
        errorbar(squeeze(mean(offer_mat(Gsub,:,k_cond,kSoc),1)),squeeze(std(offer_mat(Gsub,:,k_cond,kSoc),0,1))./sqrt(nGsub),'o',...
            'Color',colmat(k_cond,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k_cond,:))
    end
    set(gca,'XLim',[0 ntr+1],...
        'YLim',[4 12])
    xlabel('Trial')
    ylabel('Offer')
    %  legend('Cond1','Cond2','Cond3')
    title(strcat(['Social = ',num2str(kSoc-1)]));
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


SocPP = squeeze(LEARN_parametersLPP(:,:,2,2));
NoSocPP = squeeze(LEARN_parametersLPP(:,:,2,1));

[H P CI STATS] = ttest(SocPP,NoSocPP)