%% for quick and dirty plot editting

clear all
close all


load('/Users/michaelgiffin/Dropbox/RL/data_analysis_ML/BEHAVcohortwithout_priors_integr/Learning_All_2018_03_05')


%% find paths
cur_dir = '/Users/michaelgiffin/Dropbox/';
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
% [postBMC,outBMC]=VBA_groupBMC(-LEARN_LPP');



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
        
        [PE,eO,a_t, EV] = learning_models_outputs(LEARN_parametersLPP(k_sub,:,Nmodel),O,D,Nmodel);
        
        eO_mat(k_sub,:,:) = eO;
        logA_mat(k_sub,:) = a_t(end,:);
        
        EV_mat(k_sub,:,:) = EV;
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

%% Plot learning curves for offers
% h1 = figure('Units', 'pixels', ...
%     'Position', [400 200 800 250]);
% set(gcf,'Color',[1,1,1])

yL = [0 50;...
    -14 6;...
    0 5;...
    0 1;...
    0 1];


figure;
set(gcf,'Color',[1,1,1])

soccell = {'NonSocial', 'Social', };
impcell = {'Implicit', 'Explicit'};
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
        title([soccell{kSoc}, ' ', impcell{kImpexp}]);
    end
end
%


%% Plot learning curves for rewards
% h1 = figure('Units', 'pixels', ...
%     'Position', [400 200 800 250]);
% set(gcf,'Color',[1,1,1])

yL = [0 50;...
    -14 6;...
    0 5;...
    0 1;...
    0 1];

ste = @(x) (std(x))/(sqrt(length(x)));

figure;
set(gcf,'Color',[1,1,1])

soccell = {'NonSocial', 'Social', };
impcell = {'Implicit', 'Explicit'};

% impcount = 0;
% for kSoc = 1:2
%     impcount = impcount +1;
%     
%     oppcount = 0;
%     

 close all
for kImpexp = 1:2
    oppcount = 0;
    subplot(2,1,kImpexp)
    
    Gsub = impexp == (kImpexp -1);
    nGsub = sum(double(Gsub));
    
%     for kSoc = 1:2
        
        hold on
        for k_cond = 1:3
            kSoc = 1;
            
            SocCond =(kSoc-1)*3+k_cond;
            
            oppcount = oppcount+1;
            
            %             Ymod = squeeze(EV_mat(Gsub,:,SocCond));
            %             fill([(1:ntr),fliplr(1:ntr)],[mean(Ymod)+(std(Ymod)./sqrt(nGsub)),fliplr(mean(Ymod)-(std(Ymod)./sqrt(nGsub)))],...
            %                 colmat(k_cond,:),'LineStyle','none')
            %             alpha(.15)
            %             %        plot(mean(Ymod),'-',...
            %             %    'Color',colmat(k_cond,:))
            
            
            %             for kk = 1:size(beh, 1)
            %                 count    = count + 1;
            %                 colcount = colcount + 1;
            %                 bar(count, beh(kk, 1, ii)*const, mycolor{mod(colcount, 2) + 1})
            %                 h        = errorbar(count, beh(kk, 1, ii)*const, beh(kk, 2, ii)*const, 'k.', 'linewidth', 1, 'CapSize', 10);
            %                 set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % so errorbars won't be in the legend
            %             end
            %             plot([count-.5, count+.5], [nash(colcount)*const, nash(colcount)*const], '-k', 'linewidth', 5); % needs to be here to get the legend in the right place
            %             count = count+1;
            
            bar(oppcount, mean(mean(reward_mat(Gsub,:,SocCond),1)),'FaceColor', colmat(k_cond,:))
            
            h = errorbar(oppcount, mean(mean(reward_mat(Gsub,:,SocCond),1)),mean(ste(reward_mat(Gsub,:,SocCond))),symb_IE{kImpexp},...
                'Color',colmat(k_cond,:),...
                'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
                'MarkerEdgeColor',colmat(k_cond,:));
            
            kSoc = 2;
            
            oppcount = oppcount+1;
            
            SocCond =(kSoc-1)*3+k_cond;
            
            bar(oppcount, mean(mean(reward_mat(Gsub,:,SocCond),1)),'FaceColor', colmat(k_cond,:))
            
            h = errorbar(oppcount, mean(mean(reward_mat(Gsub,:,SocCond),1)),mean(ste(reward_mat(Gsub,:,SocCond))),symb_IE{kImpexp},...
                'Color',colmat(k_cond,:),...
                'MarkerFaceColor',(kSoc-1).*colmat(k_cond,:) + (2-kSoc).*[1,1,1],...
                'MarkerEdgeColor',colmat(k_cond,:));
            
            
        end
        set(gca,'XLim',[0 7],...
            'YLim',[5 16])
        xlabel('Opponent')
        ylabel('Reward')
        %  legend('Cond1','Cond2','Cond3')
        title([soccell{kSoc}, ' ', impcell{kImpexp}]);
 end
% end





%% Plot rewards in same style as paramter plots
% will have 3 subplots, one for each opponent, have the implicit condition connected with dotted
% lines and the explicit with solid lines
% opp_names = {'Opponent 1', 'Opponent 2', 'Opponent 3'};
% close all
% for k_opp = 1:3
%     
%     
%     subplot(1,3,k_opp)
%     
%     hold on
%     
%     for k_sub = 1:nImp
%         plot([.8 2.2],[mean(reward_mat(sImp(k_sub),:,k_opp+3)) mean(reward_mat(sImp(k_sub),:,k_opp))],...
%             '-.k','LineWidth',.5,'Color',.7*[1,1,1])
%     end
%         
%     plot(.8*ones(nImp,1),mean(squeeze(reward_mat(Imp,:,k_opp+3))'),'sq',...
%         'LineStyle','none',...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',colmat(k_opp,:))
%     
%     plot(2.2*ones(nImp,1),mean(reward_mat(Imp,:,k_opp)'),'sq',...
%         'LineStyle','none',...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',[1,1,1])
%     
%     errorbar([1.25 1.75],...
%         [mean(mean(reward_mat(Imp,:,k_opp+3)')),mean(mean(reward_mat(Imp,:,k_opp)'))],...
%         [mean(ste(reward_mat(Imp,:,k_opp+3)')), mean(ste(reward_mat(Imp,:,k_opp)'))],...
%         '-k','LineWidth',1,'Color',colmat(k_opp,:))
% %         [squeeze(std(LEARN_parametersLPP(Imp,5+k_param,Nmodel),0,1))./sqrt(nImp) squeeze(std(LEARN_parametersLPP(Imp,k_param,Nmodel),0,1))./sqrt(nImp)],...
% %         '-k','LineWidth',1,'Color',.2*[1,1,1])
% 
% 
%     plot(1.25,mean(mean(reward_mat(Imp,:,k_opp+3))),'sq',...
%         'LineStyle','none',...
%         'LineWidth',1,...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',colmat(k_opp,:))
%     plot(1.75,mean(mean(reward_mat(Imp,:,k_opp))),'sq',...
%         'LineStyle','none',...
%         'LineWidth',1,...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',[1,1,1])
%     
%     
%     for k_sub = 1:nExp
%         plot(2+[.8 2.2],[mean(reward_mat(sExp(k_sub),:,k_opp+3)) mean(reward_mat(sExp(k_sub),:,k_opp))],...
%             '-k','LineWidth',.5,'Color',.7*[1,1,1])
%     end
% 
%     plot(2+.8*ones(nExp,1),mean(squeeze(reward_mat(Exp,:,k_opp+3))'),'o',...
%         'LineStyle','none',...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',colmat(k_opp,:))
%     
%     plot(2+2.2*ones(nExp,1),mean(reward_mat(Exp,:,k_opp)'),'o',...
%         'LineStyle','none',...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',[1,1,1])
%     
%     errorbar(2+[1.25 1.75],...
%         [mean(mean(reward_mat(Exp,:,k_opp+3)')),mean(mean(reward_mat(Exp,:,k_opp)'))],...
%         [mean(ste(reward_mat(Exp,:,k_opp+3)')), mean(ste(reward_mat(Exp,:,k_opp)'))],...
%         '-k','LineWidth',1,'Color',colmat(k_opp,:))
% %         [squeeze(std(LEARN_parametersLPP(Exp,5+k_param,Nmodel),0,1))./sqrt(nExp) squeeze(std(LEARN_parametersLPP(Exp,k_param,Nmodel),0,1))./sqrt(nExp)],...
% %         '-k','LineWidth',1,'Color',.2*[1,1,1])
% 
% 
%     plot(2+1.25,mean(mean(reward_mat(Exp,:,k_opp+3))),'o',...
%         'LineStyle','none',...
%         'LineWidth',1,...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',colmat(k_opp,:))
%     plot(2+1.75,mean(mean(reward_mat(Exp,:,k_opp))),'o',...
%         'LineStyle','none',...
%         'LineWidth',1,...
%         'MarkerEdgeColor',colmat(k_opp,:),...
%         'MarkerFaceColor',[1,1,1])
%     
% % 
% %     hT = title(Pnames{k_param});
% %     
%     set(gca,'XLim',[0.25 4.75],...
%         'XTick', [1 2 3 4],...
%         'XTickLabel',{'S','nS'})
%     
% %     set(gca,'XLim',[0.25 4.75],...
% %         'YLim',yL(k_opp,:),...
% %         'XTick', [1 2 3 4],...
% %         'XTickLabel',{'S','nS'})
% end
% 

%% attempt at interaction plots
socpay    = [7.397778, 10.584861, 12.936944];
nonsocpay = [7.142361, 10.765694, 13.471528];
close all

figure 
for k_opp = 1:3
    hold on
      errorbar([1.25 1.75],...
        [mean(mean(reward_mat(:,:,k_opp+3)')),mean(mean(reward_mat(:,:,k_opp)'))],...
        [mean(ste(reward_mat(:,:,k_opp+3)')), mean(ste(reward_mat(:,:,k_opp)'))],...
        '-k','LineWidth',1,'Color',colmat(k_opp,:))
    
     plot(1.25,mean(mean(reward_mat(:,:,k_opp+3))),'d',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',colmat(k_opp,:),...
        'MarkerFaceColor',colmat(k_opp,:))
    plot(1.75,mean(mean(reward_mat(:,:,k_opp))),'d',...
        'LineStyle','none',...
        'LineWidth',1,...
        'MarkerEdgeColor',colmat(k_opp,:),...
        'MarkerFaceColor',[1,1,1])
    
end
ylabel('Payoff');
set(gca,'XLim',[1 2],...
    'XTick', [1.25 1.75],...
    'YLim', [6 15],...
    'XTickLabel',{'Social','NonSocial'})