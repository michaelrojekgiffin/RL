clear
close all force
clc

%% locate and Get files
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'fMRI_stimulus','UG','data');
subjects = [6:14];
nsub = numel(subjects);
ntr = 72;

%% Loop Subjects
for k = 1:nsub;
    
    
    full_data = NaN(2*ntr,7);
    sub_num = sprintf('%03.0f',(subjects(k)));
    fldir = strcat('sub',sub_num,'_prob_1.mat');
    
    flnm = fullfile(data_dir,fldir);
    load(flnm)
    
    for kcond = 1:3
        
        for kSoc = 1:2
            
            X = sub_data((sub_data(:,4) == kcond-1) & (sub_data(:,3) == kSoc-1),1:2);
            XX = sortrows(X,1);
            
            full_mat(k,kcond,:,kSoc) = XX(:,2);
        end
    end
    
end

figure;
set(gcf,'Color',[1,1,1])
colmat = [1,0,0;
    .5,0,.5;
    0,0,1];


for kSoc = 1:2
    
    subplot(1,2,kSoc)
    hold on
    
    for k = 1:3
        
        mtp = squeeze(mean(full_mat(:,k,:,kSoc),1));
        stp =  squeeze(std(full_mat(:,k,:,kSoc),0,1))./nsub;
        
        errorbar(mtp,stp,'-o',...
            'Color',colmat(k,:),...
            'MarkerFaceColor',colmat(k,:),...
            'MarkerEdgeColor',colmat(k,:));
        
    end
    
    set(gca,'YLim',[0 100],...
        'XLim',[0 20])
    xlabel('offer')
    ylabel('estimated p(accept)')
end

