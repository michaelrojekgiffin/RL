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

options     = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);

%% pre allocate
good_sub = NaN(nsub,1);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:20;
endow  = 20*ones(1,numel(offers));% parameters of the simulation

%% load priors

ntr = 25;

for k_sub = 1:nsub;
    
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
            soc_mat(k_sub,:,:,kSoc) = (kSoc-1)*ones(ntr,3);
            trial_mat(k_sub,:,:,kSoc) = repmat((1:ntr)',1,3);
            
            
            O = squeeze(offer_mat(k_sub,:,:,kSoc));
            D = squeeze(accept_mat(k_sub,:,:,kSoc));
            
            for nmodel =1:4
                n_rep           = 10;
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
end
Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

save('Learning_All_2017_12_08')