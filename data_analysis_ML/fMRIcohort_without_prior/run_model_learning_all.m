clear
close all force
clc

%% find paths
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'fMRI_stimulus','UG','data');
subjects = [3:14];
nsub = numel(subjects);
ntr = 72;

options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);

%% pre allocate
good_sub = NaN(nsub,1);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:20;
endow  = 20*ones(1,numel(offers));% parameters of the simulation

%% lrun loop

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
                    offer_mat(k_sub,:,1:3,kSoc) = [cond1(:,1),cond2(:,1),cond3(:,1)];
                    accept_mat(k_sub,:,1:3,kSoc) = [cond1(:,2),cond2(:,2),cond3(:,2)];
                    reward_mat(k_sub,:,1:3,kSoc) = [cond1(:,3),cond2(:,3),cond3(:,3)];
                    cond_mat(k_sub,:,1:3,kSoc) = [ones(24,1),2*ones(24,1),3*ones(24,1)];
                    sub_mat(k_sub,:,1:3,kSoc) = k_sub*ones(24,3);
                    soc_mat(k_sub,:,1:3,kSoc) = (kSoc-1)*ones(24,3);
                    trial_mat(k_sub,:,1:3,kSoc) = repmat((1:24)',1,3);
                case 2
                    offer_mat(k_sub,:,4:6,kSoc) = [cond1(:,1),cond2(:,1),cond3(:,1)];
                    accept_mat(k_sub,:,4:6,kSoc) = [cond1(:,2),cond2(:,2),cond3(:,2)];
                    reward_mat(k_sub,:,4:6,kSoc) = [cond1(:,3),cond2(:,3),cond3(:,3)];
                    cond_mat(k_sub,:,4:6,kSoc) = 3+[ones(24,1),2*ones(24,1),3*ones(24,1)];
                    sub_mat(k_sub,:,4:6,kSoc) = k_sub*ones(24,3);
                    soc_mat(k_sub,:,4:6,kSoc) = (kSoc-1)*ones(24,3);
                    trial_mat(k_sub,:,4:6,kSoc) = repmat((1:24)',1,3);
                    
            end
            
        end
        
        O = squeeze(offer_mat(k_sub,:,:,kSoc));
        D = squeeze(accept_mat(k_sub,:,:,kSoc));
        
        for nmodel =1:4
            n_rep           = 5;
            parameters_rep  = NaN(n_rep,5);     parametersLPP_rep  = NaN(n_rep,5);
            ll_rep          = NaN(n_rep,1);     LPP_rep          = NaN(n_rep,1);
            
            lb = [0 -15 0 0 0];          LB = [0 -Inf 0 0 0];
            ub = [5 0 5 1 1];         UB = [Inf Inf Inf 1 1];
            ddb = ub - lb;
            
            
            for k_rep = 1:n_rep
                x0 = lb + rand(1,5).*ddb;
                %standard estimation
                [parameters_rep(k_rep,1:5),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,nmodel),x0,[],[],[],[],LB,UB,[],options);
                %lalace estimation
                [parametersLPP_rep(k_rep,1:5),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,nmodel),x0,[],[],[],[],LB,UB,[],options);
            end
            
            [~,pos] = min(ll_rep);
            LEARN_parameters(k_sub,:,nmodel,kSoc)    =   parameters_rep(pos(1),:);
            LEARN_ll(k_sub,nmodel,kSoc)              =   ll_rep(pos(1),:);
            
            [~,posLPP] = min(LPP_rep);
            LEARN_parametersLPP(k_sub,:,nmodel,kSoc)      =   parametersLPP_rep(posLPP(1),:);
            LEARN_LPP(k_sub,nmodel,kSoc)                  =   LPP_rep(posLPP(1),:);
        end
        
        
    end
end
Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

save('Learning_All_2018_02_15')