clear
close all force
clc

%load('Priors_test1_2017_12_01')
load('Priors_All_2018_01_08')
%% set dir
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox

%% localize data
findnm = strfind(pwd,project_name);
%data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','rough_draft','explicit','processed');
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');

fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));

%% load and set params

offers  = 0:1:20;

good_sub = NaN(nsub,1);

%% subject loop
for k_sub = 1:nsub;
    
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    flnmDict = fullfile(data_dir,strcat('Sub_',sub_nm,'_Dictator.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        good_sub(k_sub) = 1;
        
        for kSoc = 1:3
            %% Get UG data
            
            switch kSoc
                case 1
                    load(flnmNonSoc)
                case 2
                    load(flnmSoc)
                case 3
                    load(flnmDict)
            end
            
            %% get the distribs
            PA_sub(k_sub,:,kSoc)     = logitp([squeeze(parametersLPP(k_sub,1,kSoc)),squeeze(parametersLPP(k_sub,2,kSoc))],offers);            % compute proba of accepting the offers given current model
            EV_sub(k_sub,:,kSoc)    = (endow - offers).* squeeze(PA_sub(k_sub,:,kSoc)) ;
            
            for k = 1:length(subdata)
                pc(k)     = exp(parametersLPP(k_sub,3,kSoc).*EV_sub(k_sub,subdata(k,1)+1,kSoc)) ./ sum(exp(parametersLPP(k_sub,3,kSoc).*EV_sub(k_sub,subdata(k,1:2)+1,kSoc)));                                % resample Offer in pdf (="soft-max")
            end
            
            
            ch  = subdata(:,3) == subdata(:,1);
            chO = subdata(:,1:2);
            
            for kk = 1:length(offers)
                behav_sub(k_sub,kk,kSoc) = nanmean([mean(ch(chO(:,1)==offers(kk))),1-mean(ch(chO(:,2)==offers(kk)))],2);
                pc_sub(k_sub,kk,kSoc) = nanmean([mean(pc(chO(:,1)==offers(kk))),1-mean(pc(chO(:,2)==offers(kk)))],2);
            end
            
            
        end
    end
    
end

Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

figure;
set(gcf,'Color',[1,1,1])

subplot(2,2,1)
hold on
errorbar(offers,squeeze(mean(behav_sub(Gsub,:,3))),squeeze(std(behav_sub(Gsub,:,3)))./sqrt(nGsub),'o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',0.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(behav_sub(Gsub,:,1))),squeeze(std(behav_sub(Gsub,:,1)))./sqrt(nGsub),'o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(behav_sub(Gsub,:,2))),squeeze(std(behav_sub(Gsub,:,2)))./sqrt(nGsub),'o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',.5.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])

plot(offers,squeeze(mean(pc_sub(Gsub,:,3))),'--',...
    'Color',[0,0,0])
plot(offers,squeeze(mean(pc_sub(Gsub,:,1))),'-',...
    'Color',.7.*[1,1,1])
plot(offers,squeeze(mean(pc_sub(Gsub,:,2))),'--',...
    'Color',.3*[0,0,0])

set(gca,'XLim',[0 21],...
    'YLim',[0 1])
xlabel('offers')
ylabel('% Chosen')
legend('Dict','NonSoc','Soc')


subplot(2,2,2)
hold on

mtp = squeeze(mean(parametersLPP,1));
stp = squeeze(std(parametersLPP,0,1))./sqrt(nGsub);
bar(mtp)
legend('NonSoc','Soc','Dict')

subplot(2,2,3)
hold on
errorbar(offers,squeeze(mean(EV_sub(Gsub,:,3))),squeeze(std(EV_sub(Gsub,:,3)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',0.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(EV_sub(Gsub,:,1))),squeeze(std(EV_sub(Gsub,:,1)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(EV_sub(Gsub,:,2))),squeeze(std(EV_sub(Gsub,:,2)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',.5.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])

set(gca,'XLim',[0 21],...
    'YLim',[0 8])
xlabel('offers')
ylabel('Estimated Expected value')
legend('Dict','NonSoc','Soc')


subplot(2,2,4)
hold on
errorbar(offers,squeeze(mean(PA_sub(Gsub,:,3))),squeeze(std(PA_sub(Gsub,:,3)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',0.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(PA_sub(Gsub,:,1))),squeeze(std(PA_sub(Gsub,:,1)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(PA_sub(Gsub,:,2))),squeeze(std(PA_sub(Gsub,:,2)))./sqrt(nGsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',.5.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[0 21],...
    'YLim',[0 1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')
legend('Dict','NonSoc','Soc')







figure
set(gcf,'Color',[1,1,1])
for k = 1:3
    subplot(2,3,k)
    plot(parameters(Gsub,k,1),parametersLPP(Gsub,k,1),'o',...
        'MarkerFaceColor',.5*[1,1,1],...
        'MarkerEdgeColor',[0,0,0])
    xlabel('LogLik')
    ylabel('Laplace LPP')
    subplot(2,3,3+k)
    plot(parameters(Gsub,k,2),parametersLPP(Gsub,k,2),'o',...
        'MarkerFaceColor',.9*[1,1,1],...
        'MarkerEdgeColor',[0,0,0])
    xlabel('LogLik')
    ylabel('Laplace LPP')
end

[H,P,CI,STATS] = ttest(squeeze(parametersLPP(:,:,1)),squeeze(parametersLPP(:,:,2)));

