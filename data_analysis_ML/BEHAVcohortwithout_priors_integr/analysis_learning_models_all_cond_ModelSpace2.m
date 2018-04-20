clear
close all force
clc


%% load priors
load('Learning_All_2018_03_05')

%% find paths
cur_dir = pwd;
% project_name = 'RL_PreyPredator';
project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');
fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));
nsub = length(fl_list);



%% some params
nfpm = [5,...
    6,6,6,6,6,...
    7,7,7,7,7,7,7,7,7,7,...
    8,8,8,8,8,8,8,8,8,8,...
    9,9,9,9,9,...
    10];
Nmodel = 32;
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
X = squeeze(LEARN_LPP(:,:));
for k_mod= 1:17
    bic(:,k_mod)=-2*-X(:,k_mod) + nfpm(k_mod)*log(2*3*25); % l2 is already positive
end
% [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
[postBMC,outBMC]=VBA_groupBMC(-LEARN_LPP');


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
        
        
        SOC = subdata(subdata(:,4) == 0,:);
        cond1 = SOC(SOC(:,5)==1,:);
        cond2 = SOC(SOC(:,5)==2,:);
        cond3 = SOC(SOC(:,5)==3,:);
        offer_mat(k_sub,:,1:3) = [cond1(1:ntr,1),cond2(1:ntr,1),cond3(1:ntr,1)];
        accept_mat(k_sub,:,1:3) = [cond1(1:ntr,2),cond2(1:ntr,2),cond3(1:ntr,2)];
        reward_mat(k_sub,:,1:3) = [cond1(1:ntr,3),cond2(1:ntr,3),cond3(1:ntr,3)];
        cond_mat(k_sub,:,1:3) = [ones(ntr,1),2*ones(ntr,1),3*ones(ntr,1)];
        sub_mat(k_sub,:,1:3) = k_sub*ones(ntr,3);
        soc_mat(k_sub,:,1:3) = zeros(ntr,3);
        trial_mat(k_sub,:,1:3) = repmat((1:ntr)',1,3);
        
        SOC = subdata(subdata(:,4) == 1,:);
        cond1 = SOC(SOC(:,5)==1,:);
        cond2 = SOC(SOC(:,5)==2,:);
        cond3 = SOC(SOC(:,5)==3,:);
        offer_mat(k_sub,:,4:6) = [cond1(1:ntr,1),cond2(1:ntr,1),cond3(1:ntr,1)];
        accept_mat(k_sub,:,4:6) = [cond1(1:ntr,2),cond2(1:ntr,2),cond3(1:ntr,2)];
        reward_mat(k_sub,:,4:6) = [cond1(1:ntr,3),cond2(1:ntr,3),cond3(1:ntr,3)];
        cond_mat(k_sub,:,4:6) = [ones(ntr,1),2*ones(ntr,1),3*ones(ntr,1)];
        sub_mat(k_sub,:,4:6) = k_sub*ones(ntr,3);
        soc_mat(k_sub,:,4:6) = ones(ntr,3);
        trial_mat(k_sub,:,4:6) = repmat((1:ntr)',1,3);
        
        
        
        O = squeeze(offer_mat(k_sub,:,:));
        D = squeeze(accept_mat(k_sub,:,:));
        
        [PE,eO,a_t] = learning_models_outputs(LEARN_parametersLPP(k_sub,:,Nmodel),O,D,Nmodel);
        
        eO_mat(k_sub,:,:) = eO;
        logA_mat(k_sub,:) = a_t(end,:);
        %         logBmat(k_sub,:) = repmat(LEARN_parametersLPP(k_sub,3,Nmodel,kSoc),1,3);
        %
        %         for kcond = 1:3
        %             PA_sub(k_sub,:,kcond,kSoc)     = logitp([squeeze(logA_mat(k_sub,kcond,kSoc))',squeeze(logBmat(k_sub,kcond,kSoc))'],offers);            % compute proba of accepting the offers given current model
        %             EV_sub(k_sub,:,kcond,kSoc)    = (endow - offers).* squeeze(PA_sub(k_sub,:,kcond,kSoc)) ;
        %         end
    end
end





colmat = [1,0,0;.5,0,.5;0,0,1];
symb_IE = {'sq','o','d'};
Imp = impexp == 0;
sImp = find(impexp == 0);
nImp = sum(double(Imp));
Exp = impexp == 1;
sExp = find(impexp == 1);
nExp =  sum(double(Exp));

%% Plot learning curves
h1 = figure('Units', 'pixels', ...
    'Position', [400 200 800 250]);
set(gcf,'Color',[1,1,1])

yL = [0 50;...
    -14 6;...
    0 5;...
    0 1;...
    0 1];

Pnames = {'\beta','a_0','b_0','\alpha_1','\alpha_2'};

%% Plot Average Params
for k_param = 1:5
    
    
    subplot(1,5,k_param)
    hold on
    for k_sub = 1:nImp
        plot([.8 2.2],[LEARN_parametersLPP(sImp(k_sub),5+k_param,Nmodel) LEARN_parametersLPP(sImp(k_sub),k_param,Nmodel)],...
            '-k','LineWidth',.5,'Color',.7*[1,1,1])
    end
    plot(.8*ones(nImp,1),squeeze(LEARN_parametersLPP(Imp,5+k_param,Nmodel)),'sq',...
        'LineStyle','none',...
        'MarkerEdgeColor',.5.*[1,1,1],...
        'MarkerFaceColor',.6.*[1,1,1])
    plot(2.2*ones(nImp,1),squeeze(LEARN_parametersLPP(Imp,k_param,Nmodel)),'sq',...
        'LineStyle','none',...
        'MarkerEdgeColor',.7.*[1,1,1],...
        'MarkerFaceColor',[1,1,1])
    
    errorbar([1.25 1.75],...
        [squeeze(mean(LEARN_parametersLPP(Imp,5+k_param,Nmodel),1)) squeeze(mean(LEARN_parametersLPP(Imp,k_param,Nmodel),1))],...
        [squeeze(std(LEARN_parametersLPP(Imp,5+k_param,Nmodel),0,1))./sqrt(nImp) squeeze(std(LEARN_parametersLPP(Imp,k_param,Nmodel),0,1))./sqrt(nImp)],...
        '-k','LineWidth',1,'Color',.2*[1,1,1])
    plot(1.25,squeeze(mean(LEARN_parametersLPP(Imp,5+k_param,Nmodel),1)),'sq',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',.0.*[1,1,1])
    plot(1.75,squeeze(mean(LEARN_parametersLPP(Imp,k_param,Nmodel),1)),'sq',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',[1,1,1])
    
    
    
    for k_sub = 1:nExp
        plot(2+[.8 2.2],[LEARN_parametersLPP(sExp(k_sub),5+k_param,Nmodel) LEARN_parametersLPP(sExp(k_sub),k_param,Nmodel)],...
            '-k','LineWidth',.5,'Color',.7*[1,1,1])
    end
    plot(2+.8*ones(nExp,1),squeeze(LEARN_parametersLPP(Exp,5+k_param,Nmodel)),'o',...
        'LineStyle','none',...
        'MarkerEdgeColor',.5.*[1,1,1],...
        'MarkerFaceColor',.6.*[1,1,1])
    plot(2+2.2*ones(nExp,1),squeeze(LEARN_parametersLPP(Exp,k_param,Nmodel)),'o',...
        'LineStyle','none',...
        'MarkerEdgeColor',.7.*[1,1,1],...
        'MarkerFaceColor',[1,1,1])
    
    errorbar(2+[1.25 1.75],...
        [squeeze(mean(LEARN_parametersLPP(Exp,5+k_param,Nmodel),1)) squeeze(mean(LEARN_parametersLPP(Exp,k_param,Nmodel),1))],...
        [squeeze(std(LEARN_parametersLPP(Exp,5+k_param,Nmodel),0,1))./sqrt(nExp) squeeze(std(LEARN_parametersLPP(Exp,k_param,Nmodel),0,1))./sqrt(nExp)],...
        '-k','LineWidth',1,'Color',.2*[1,1,1])
    plot(2+1.25,squeeze(mean(LEARN_parametersLPP(Exp,5+k_param,Nmodel),1)),'o',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',.0.*[1,1,1])
    plot(2+1.75,squeeze(mean(LEARN_parametersLPP(Exp,k_param,Nmodel),1)),'o',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',[0,0,0],...
        'MarkerFaceColor',[1,1,1])
    
    hT = title(Pnames{k_param});
    
    set(gca,'XLim',[0.25 4.75],...
        'YLim',yL(k_param,:),...
        'XTick', [1 2 3 4],...
        'XTickLabel',{'S','nS'})
end

%
 pNSoc = squeeze(LEARN_parametersLPP(:,1:5,Nmodel));
 pSoc = squeeze(LEARN_parametersLPP(:,6:10,Nmodel));
 [H P CI STATS] = ttest(pNSoc,pSoc)
 
% pNSoc = [easy_sub{1,1};easy_sub{1,2}];
% pSoc  = [easy_sub{2,1};easy_sub{2,2}];
% 
%
% pImp = [easy_sub{1,1}+easy_sub{2,1}]./2;
% pExp  = [easy_sub{1,2}+easy_sub{2,2}]./2;
% [H P CI STATS] = ttest2(pImp,pExp)
%
% save(['easy_sub' date], 'easy_sub');
% % % SocExp vs. NonSocExp Beta
% % [H P CI STATS] = ttest2(easy_sub{1,2}(:,1),easy_sub{2,2}(:,1))
% % % intercept
% % [H P CI STATS] = ttest2(easy_sub{1,2}(:,2),easy_sub{2,2}(:,2))
% % %lr 1
% % [H P CI STATS] = ttest2(easy_sub{1,2}(:,4),easy_sub{2,2}(:,4))

%% Plot learning curves

figure;
set(gcf,'Color',[1,1,1])


for kSoc = 1:2
    
    for kImpexp = 1:2
        
        subplot(2,2,(kImpexp-1)*2 + kSoc)
        Gsub = impexp == (kImpexp -1);
        nGsub = sum(double(Gsub));
        
        hold on
        for k_cond = 1:3
            SocCond =(kSoc-1)*3+k_cond;
            Ymod = squeeze(eO_mat(Gsub,:,SocCond));
            fill([(1:ntr),fliplr(1:ntr)],[mean(Ymod)+(std(Ymod)./sqrt(nGsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nGsub)))],...
                colmat(k_cond,:),'LineStyle','none')
            alpha(.15)
            %        plot(mean(Ymod),'-',...
            %    'Color',colmat(k_cond,:))
            
            errorbar(squeeze(mean(offer_mat(Gsub,:,SocCond),1)),squeeze(std(offer_mat(Gsub,:,SocCond),0,1))./sqrt(nGsub),symb_IE{kImpexp},...
                'Color',colmat(k_cond,:),...
                'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
                'MarkerEdgeColor',colmat(k_cond,:))
        end
        set(gca,'XLim',[0 ntr+1],...
            'YLim',[3.5 9.5])
        xlabel('Trial')
        ylabel('Offer')
        %  legend('Cond1','Cond2','Cond3')
        title(strcat(['Social = ',num2str(kSoc-1) '; ImpExp == ',num2str(kImpexp-1)]));
    end
end
%
% %% Plot observed acceptance rate
% figure
% set(gcf,'Color',[1,1,1])
%
% for kSoc = 1:2
%
%     subplot(1,2,kSoc)
%     hold on
%
%     for k_cond = 1:3
%
%         offer = squeeze(offer_mat(Gsub,:,k_cond,kSoc));
%         offer = offer(:);
%         accept = squeeze(accept_mat(Gsub,:,k_cond,kSoc));
%         accept = accept(:);
%
%         nOFF = unique(offer);
%         off_rate =[];
%         for koff = 1:length(nOFF);
%             off_rate(koff) = mean(accept(offer == koff));
%         end
%
%         % plot final model acceptance estimations
%         Ymod = squeeze(PA_sub(Gsub,:,k_cond,kSoc));
%         fill([(offers),fliplr(offers)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nsub)))],...
%             colmat(k_cond,:),'LineStyle','none')
%         alpha(.15)
%         plot(offers,mean(Ymod),'-',...
%             'Color',colmat(k_cond,:))
%
%         % plot observed acceptance frequency (~estimation if subjects are optimal)
%         plot(nOFF,off_rate,'o',...
%             'Color',colmat(k_cond,:),...
%             'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
%             'MarkerEdgeColor',colmat(k_cond,:))
%
%         set(gca,'YLim',[0 1])
%     end
%
%     set(gca,'XLim',[0 20],...
%         'YLim',[0 1])
%     xlabel('Offer')
%     ylabel('Acceptance rate')
%     title(strcat(['Social = ',num2str(kSoc-1)]))
%
% end
%
% %% Plot expected reward
% figure
% set(gcf,'Color',[1,1,1])
%
% for kSoc = 1:2
%
%     subplot(1,2,kSoc)
%     hold on
%
%     for k_cond = 1:3
%
%         offer = squeeze(offer_mat(Gsub,:,k_cond,kSoc));
%         offer = offer(:);
%         reward = squeeze(reward_mat(Gsub,:,k_cond,kSoc));
%         reward = reward(:);
%
%         nOFF = unique(offer);
%         off_rate =[];
%         for koff = 1:length(nOFF);
%             off_rate(koff) = mean(reward(offer == koff));
%         end
%
%         % plot model
%         Ymod = squeeze(EV_sub(Gsub,:,k_cond,kSoc));
%         plot(offers,mean(Ymod),'-',...
%             'Color',colmat(k_cond,:))
%         fill([(offers),fliplr(offers)],[mean(Ymod)+(std(Ymod)./sqrt(nsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nsub)))],...
%             colmat(k_cond,:),'LineStyle','none')
%         alpha(.15)
%
%         plot(nOFF,off_rate,'o',...
%             'Color',colmat(k_cond,:),...
%             'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
%             'MarkerEdgeColor',colmat(k_cond,:))
%
%
%
%         set(gca,'YLim',[0 1])
%     end
%
%     set(gca,'XLim',[0 20],...
%         'YLim',[0 20])
%     xlabel('Offer')
%     ylabel('Expected value')
%     title(strcat(['Social = ',num2str(kSoc-1)]))
%
% end
%
%
%
