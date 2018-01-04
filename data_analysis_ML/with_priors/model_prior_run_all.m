clear
close all force
clc



%% set dir
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox

%% set optioms
rng('shuffle')
options         = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);


%% localize data
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed');

fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));

%% useful functions
offers  = 0:1:20;
endow   = 20*ones(1,numel(offers));
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


nsub = length(fl_list);

%% subject loop
gsub =0;
for k_sub = 1:nsub;
    
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        gsub = gsub+1;
        prior_sub_names{gsub} = sub_nm;
        
        for kSoc = 1:2
            %% Get UG data
            
            switch kSoc
                case 1
               load(flnmNonSoc)
                case 2
               load(flnmSoc)
            end
            
            ch  = subdata(:,3) == subdata(:,1);
            chO = subdata(:,1:2);
            
            for kk = 1:20
               behav(k_sub,kk) = mean(ch(chO(:,1)==kk)) + 1-mean(ch(chO(:,2)==kk));          
            end
            
            full_o(k_sub,:,:,kSoc) = offers;
            full_ch(k_sub,:,kSoc)  = ch;
            
            %% fit the prior
            n_rep           = 5;
            parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
            ll_rep          = NaN(n_rep,1);     LPP_rep             = NaN(n_rep,1);
            
            lb = [-20 0 0];     LB = [-Inf 0 0];
            ub = [10 20 20];    UB = [Inf Inf Inf];
            
            for k_rep = 1:n_rep
                x0 = [-5*rand() 5*rand()  5*rand()];
                
                % standard estim
                [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) estimate_priors_bin(x,chO,ch),x0,[],[],[],[],lb,ub,[],options);
                % laplace approximation
                [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_priors(x,chO,ch),x0,[],[],[],[],LB,UB,[],options);
            end
            [~,pos]                  = min(ll_rep);
            parameters(k_sub,:,kSoc) = parameters_rep(pos(1),:);
            ll(k_sub)                = ll_rep(pos(1),:);
            
            [~,posLPP]                  = min(LPP_rep);
            parametersLPP(k_sub,:,kSoc) = parametersLPP_rep(posLPP(1),:);
            LPP(k_sub)                  = LPP_rep(posLPP(1),:);
            
            %% get the distribs
            PA_sub(k_sub,:,kSoc)     = logitp([squeeze(parameters(k_sub,1,kSoc)),squeeze(parameters(k_sub,2,kSoc))],offers);            % compute proba of accepting the offers given current model
            EV_sub(k_sub,:,kSoc)    = (endow - offers).* squeeze(PA_sub(k_sub,:,kSoc)) ;                                   % compute EV of the offers given current model
        end
    end
end


save('Priors_All_2017_12_07')


figure;
set(gcf,'Color',[1,1,1])


subplot(1,2,1)
hold on
errorbar(offers,squeeze(mean(EV_sub(:,:,1))),squeeze(std(EV_sub(:,:,1)))./sqrt(nsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(EV_sub(:,:,2))),squeeze(std(EV_sub(:,:,2)))./sqrt(nsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',.5.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[0 21],...
    'YLim',[0 8])
xlabel('offers')
ylabel('Estimated Expected value')
legend('NonSoc','Soc')

subplot(1,2,2)
hold on
errorbar(offers,squeeze(mean(PA_sub(:,:,1))),squeeze(std(PA_sub(:,:,1)))./sqrt(nsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
errorbar(offers,squeeze(mean(PA_sub(:,:,2))),squeeze(std(PA_sub(:,:,2)))./sqrt(nsub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',.5.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[0 21],...
    'YLim',[0 1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')
legend('NonSoc','Soc')



figure
set(gcf,'Color',[1,1,1])
for k = 1:3
subplot(2,3,k)
plot(parameters(:,k,1),parametersLPP(:,k,1),'o')
xlabel('LogLik')
ylabel('Laplace LPP')
subplot(2,3,3+k)
plot(parameters(:,k,2),parametersLPP(:,k,2),'o')
xlabel('LogLik')
ylabel('Laplace LPP')
end
            
            

