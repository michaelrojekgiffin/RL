clear
close all force
clc

%% locate and Get files
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');
fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));
nsub = length(fl_list);

ntr = 25;
%% Loop Subjects
for k = 1:nsub;
    
    sub_nm = fl_list(k).name(5:12);
    
    
    %% Get UG data
    flnm = fullfile(data_dir,strcat('Sub_',sub_nm,'_UG.mat'));
    load(flnm)
    
    for kC = 1:2
        SOC = subdata(subdata(:,4) == kC-1,:);
        cond1 = SOC(SOC(:,5)==1,:);
        cond2 = SOC(SOC(:,5)==2,:);
        cond3 = SOC(SOC(:,5)==3,:);
        offer_mat(k,:,:,kC) = [cond1(1:ntr,1),cond2(1:ntr,1),cond3(1:ntr,1)];
        accept_mat(k,:,:,kC) = [cond1(1:ntr,2),cond2(1:ntr,2),cond3(1:ntr,2)];
        reward_mat(k,:,:,kC) = [cond1(1:ntr,3),cond2(1:ntr,3),cond3(1:ntr,3)];
        cond_mat(k,:,:,kC) = [ones(ntr,1),2*ones(ntr,1),3*ones(ntr,1)];
        sub_mat(k,:,:,kC) = k*ones(ntr,3);
        soc_mat(k,:,:,kC) = (kC-1)*ones(ntr,3);
        trial_mat(k,:,:,kC) = repmat((1:ntr)',1,3);
        
    end
    
end



colmat = [1,0,0;.5,0,.5;0,0,1];

%% Plot learning curves

figure;
set(gcf,'Color',[1,1,1])


for kSoc = 1:2

subplot(1,2,kSoc)
hold on
for k = 1:3
    errorbar(squeeze(mean(offer_mat(:,:,k,kSoc),1)),squeeze(std(offer_mat(:,:,k,kSoc),0,1))./sqrt(nsub),'-o',...
        'Color',colmat(k,:),...
        'MarkerFaceColor',(kSoc-1).*colmat(k,:) + (2-kSoc).*[1,1,1],...
        'MarkerEdgeColor',colmat(k,:))
end
set(gca,'XLim',[0 ntr+1],...
    'YLim',[4 10])
xlabel('Trial')
ylabel('Offer')
legend('Cond1','Cond2','Cond3')
title(strcat(['Social = ',num2str(kSoc -1)]));
end


%% Plot Average offers


X_SOC = squeeze(mean(offer_mat(:,:,:,2),2));
X_non_SOC = squeeze(mean(offer_mat(:,:,:,1),2));

mtp_SOC = squeeze(mean(X_SOC,1));
stp_SOC = squeeze(std(X_SOC,0,1))./sqrt(nsub);

mtp_non_SOC = squeeze(mean(X_non_SOC,1));
stp_non_SOC = squeeze(std(X_non_SOC,0,1))./sqrt(nsub);

figure;
set(gcf,'Color',[1,1,1])
subplot(1,2,1)
hold on
for k = 1:3
    bar(k,mtp_non_SOC(:,k),...
        'FaceColor',[1,1,1],...
        'EdgeColor',colmat(k,:))
    errorbar(k,mtp_non_SOC(:,k),stp_non_SOC(:,k),'k','LineStyle','none')
end

set(gca,'XTick',1:3,...
    'YLim',[4 10])
xlabel('COND')
ylabel('Average Offer')
title('Social = 0')

subplot(1,2,2)
hold on


for k = 1:3
    bar(k,mtp_SOC(:,k),...
        'FaceColor',colmat(k,:),...
        'EdgeColor',[0,0,0])
    errorbar(k,mtp_SOC(:,k),stp_SOC(:,k),'k','LineStyle','none')
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
    
    for k = 1:3
        
        offer = squeeze(offer_mat(:,:,k,kSoc));
        offer = offer(:);
        accept = squeeze(accept_mat(:,:,k,kSoc));
        accept = accept(:);
        
        nOFF = unique(offer);
        off_rate =[];
        for koff = 1:length(nOFF);
            off_rate(koff) = mean(accept(offer == koff));
        end
        
        plot(nOFF,off_rate,'-o',...
            'Color',colmat(k,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k,:))
        
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
    
    for k = 1:3
        
        offer = squeeze(offer_mat(:,:,k,kSoc));
        offer = offer(:);
        reward = squeeze(reward_mat(:,:,k,kSoc));
        reward = reward(:);
        
        nOFF = unique(offer);
        off_rate =[];
        for koff = 1:length(nOFF);
            off_rate(koff) = mean(reward(offer == koff));
        end
        
        plot(nOFF,off_rate,'-o',...
            'Color',colmat(k,:),...
            'MarkerFaceColor',(kSoc-1).*colmat(k,:) + (2-kSoc).*[1,1,1],...
            'MarkerEdgeColor',colmat(k,:))
        
        set(gca,'YLim',[0 1])
    end
        
set(gca,'XLim',[0 20],...
    'YLim',[0 20])
xlabel('Offer')
ylabel('Expected value')
title(strcat(['Social = ',num2str(kSoc-1)]))

end



%% Simple ANOVA
% SUB = repmat((1:60)',1,3);
% sub = repmat(SUB,2,1);
% COND = [ones(60,1),2*ones(60,1),3*ones(60,1)];
% cond = repmat(COND,2,1);
%
% manip = [ones(60,3);zeros(60,3)];
% y = [X_SOC;X_non_SOC];
%
% vnames = {'social','cond','sub'};
%  [p,table,stats]=anovan(y(:),{manip(:),cond(:),sub(:)},'random',3,'model','interaction','varnames',vnames);
%  [p_cont,table_cont,stats_cont]=anovan(y(:),{manip(:),cond(:),sub(:)},'random',3,'continuous',2,'model','interaction','varnames',vnames);
%
%
%% FULL mixed-effect model
%  offer = offer_mat(:);
%  sub = sub_mat(:);
%  cond = cond_mat(:);
%  soc = soc_mat(:);
%  trial = trial_mat(:);
%  vnames = {'social','cond','trial','sub'};
%  [p_cont,table_cont,stats_cont]=anovan(offer,{soc,cond,trial,sub},'random',4,'continuous',[2,3],'model','full','varnames',vnames);
% 




