clear
close all force
clc

%% find paths
cur_dir         = pwd;
project_name    = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox
findnm          = strfind(pwd,project_name);
data_dir        = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');
fl_list         = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));
nsub            = length(fl_list);

options         = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);

%% pre allocate
good_sub            = NaN(nsub,1);
LEARN_parameters    = NaN(100,10,32);
LEARN_ll            = NaN(100,32);
LEARN_parametersLPP = NaN(100,10,32);
LEARN_LPP           = NaN(100,32);

%% pre-allocate
ntr         = 25;
offer_mat   = NaN(100,ntr,6);
accept_mat  = NaN(100,ntr,6);
reward_mat  = NaN(100,ntr,6);
cond_mat    = NaN(100,ntr,6);
sub_mat     = NaN(100,ntr,6);
soc_mat     = NaN(100,ntr,6);
trial_mat   = NaN(100,ntr,6);

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:20;
endow  = 20*ones(1,numel(offers));% parameters of the simulation

for k_sub = 1:nsub;
    
    %% check sub
    sub_nm      = fl_list(k_sub).name(5:12);
    flnmSoc     = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc  = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        good_sub(k_sub) = 1;
        %% Get UG data
        flnm = fullfile(data_dir,strcat('Sub_',sub_nm,'_UG.mat'));
        load(flnm)
        
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
        S = squeeze(soc_mat(k_sub,:,:));
        
        %% model fitting
        n_rep           = 10;
        for nmodel = 1:32
            
            LL_bds = [0 -15 0 0 0;...
                5 0 5 1 1];
            LPP_bds = [0 -Inf 0 0 0;...
                Inf Inf Inf 1 1];
            
            switch nmodel
                
                case 1
                    np = 5;
                    lb = LL_bds(1,:);                   LB = LPP_bds(1,:);
                    ub = LL_bds(2,:);                   UB = LPP_bds(2,:);
                case 2
                    np = 6;
                    lb = [LL_bds(1,:),LL_bds(1,1)];     LB = [LPP_bds(1,:),LPP_bds(1,1)];
                    ub = [LL_bds(2,:),LL_bds(2,1)];     UB = [LPP_bds(2,:),LPP_bds(2,1)];
                case 3
                    np = 6;
                    lb = [LL_bds(1,:),LL_bds(1,2)];     LB = [LPP_bds(1,:),LPP_bds(1,2)];
                    ub = [LL_bds(2,:),LL_bds(2,2)];     UB = [LPP_bds(2,:),LPP_bds(2,2)];
                case 4
                    np = 6;
                    lb = [LL_bds(1,:),LL_bds(1,3)];     LB = [LPP_bds(1,:),LPP_bds(1,3)];
                    ub = [LL_bds(2,:),LL_bds(2,3)];     UB = [LPP_bds(2,:),LPP_bds(2,3)];
                case 5
                    np = 6;
                    lb = [LL_bds(1,:),LL_bds(1,4)];     LB = [LPP_bds(1,:),LPP_bds(1,4)];
                    ub = [LL_bds(2,:),LL_bds(2,4)];     UB = [LPP_bds(2,:),LPP_bds(2,4)];
                case 6
                    np = 6;
                    lb = [LL_bds(1,:),LL_bds(1,5)];     LB = [LPP_bds(1,:),LPP_bds(1,5)];
                    ub = [LL_bds(2,:),LL_bds(2,5)];     UB = [LPP_bds(2,:),LPP_bds(2,5)];
                case 7
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2])];
                case 8
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[1 3])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 3])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 3])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 3])];
                case 9
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[1 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 4])];
                case 10
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[1 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 5])];
                case 11
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[2 3])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 3])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 3])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 3])];
                case 12
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[2 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 4])];
                case 13
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[2 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 5])];
                case 14
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[3 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[3 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[3 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[3 4])];
                case 15
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[3 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[3 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[3 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[3 5])];
                case 16
                    np = 7;
                    lb = [LL_bds(1,:),LL_bds(1,[4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[4 5])];
                case 17
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 3])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 3])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 3])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 3])];
                case 18
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 4])];
                case 19
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 5])];
                case 20
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 3 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 3 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 3 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 3 4])];
                case 21
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 3 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 3 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 3 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 3 5])];
                case 22
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[1 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 4 5])];
                case 23
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[2 3 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 3 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 3 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 3 4])];
                case 24
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[2 3 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 3 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 3 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 3 5])];
                case 25
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[2 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 4 5])];
                case 26
                    np = 8;
                    lb = [LL_bds(1,:),LL_bds(1,[3 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[3 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[3 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[3 4 5])];
                case 27
                    np = 9;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 3 4])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 3 4])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 3 4])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 3 4])];
                case 28
                    np = 9;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 3 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 3 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 3 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 3 5])];
                case 29
                    np = 9;
                    lb = [LL_bds(1,:),LL_bds(1,[1 2 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 2 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 2 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 2 4 5])];
                case 30
                    np = 9;
                    lb = [LL_bds(1,:),LL_bds(1,[1 3 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[1 3 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[1 3 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[1 3 4 5])];
                case 31
                    np = 9;
                    lb = [LL_bds(1,:),LL_bds(1,[2 3 4 5])];     LB = [LPP_bds(1,:),LPP_bds(1,[2 3 4 5])];
                    ub = [LL_bds(2,:),LL_bds(2,[2 3 4 5])];     UB = [LPP_bds(2,:),LPP_bds(2,[2 3 4 5])];
                case 32
                    np = 10;
                    lb = repmat(LL_bds(1,:),1,2);       LB = repmat([0 -Inf 0 0 0],1,2);
                    ub = repmat(LL_bds(2,:),1,2);       UB = repmat([Inf Inf Inf 1 1],1,2);
                    
            end
            
            ddb = ub - lb;
            parameters_rep  = NaN(n_rep,np);     parametersLPP_rep  = NaN(n_rep,np);
            ll_rep          = NaN(n_rep,1);      LPP_rep          = NaN(n_rep,1);
            
            for k_rep = 1:n_rep
                x0 = lb + rand(1,np).*ddb;
                %standard estimation               
                [parameters_rep(k_rep,1:np),ll_rep(k_rep,1)]=fmincon(@(x) learning_models_estim(x,O,D,nmodel),x0,[],[],[],[],LB,UB,[],options);
                %lalace estimation
                [parametersLPP_rep(k_rep,1:np),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning(x,O,D,nmodel),x0,[],[],[],[],LB,UB,[],options);
            end
            
            [~,pos] = min(ll_rep);
            LEARN_parameters(k_sub,1:np,nmodel)     =   parameters_rep(pos(1),:);
            LEARN_ll(k_sub,nmodel)                  =   ll_rep(pos(1),:);
            
            [~,posLPP] = min(LPP_rep);
            LEARN_parametersLPP(k_sub,1:np,nmodel)  =   parametersLPP_rep(posLPP(1),:);
            LEARN_LPP(k_sub,nmodel)                  =   LPP_rep(posLPP(1),:);
        end
        
    end
end
Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

save('Learning_All_2018_03_05')