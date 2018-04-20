% to make some model free plots corroborating the statistics I've found in
% R

clear
close all force
clc


% standard error function
ste = @(x) (std(x))/(sqrt(length(x)));
% ste = @(x) (std(x))/(sqrt(length(x)));

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% find paths and import data
% this gets all subjects into one big matrix, not super elegant but closer
% to what I'm used to in R so should help speed things along
cur_dir = pwd;
%  project_name = 'RL_PreyPredator';
project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);

data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed_R', 'UG.mat');

load(data_dir);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% for the UG datafile
% column 1: offer
% column 2: accepted (1 for yes, 0 for no)
% column 3: payment
% column 4: social (1 for social, 0 for non-social)
% column 5: opponent, 1 for starting endowment of 0, 2 for starting endowment of 10, and 3 for starting endowment of 20
% column 6: explicit (1 for explicit, 0 for implicit)
% column 7: subject name
% column 8: trial
% column 9: subject age
all_header = {'offer', 'accepted', 'payoff', 'social', 'opponent', 'explicit', 'sub_name', 'trial', 'age'};
%% format data in nice way
ntr = 25;

soc     = UG(UG(:,4) == 1, :);

soc_imp = soc(soc(:,6) == 0, :);
soc_exp = soc(soc(:,6) == 1, :);

nonsoc  = UG(UG(:,4) == 0, :);

nonsoc_exp = nonsoc(nonsoc(:,6) == 1, :);
nonsoc_imp = nonsoc(nonsoc(:,6) == 0,  :);

explicit = UG(UG(:,6) == 1, :);
implicit = UG(UG(:,6) == 0, :);

% for ugm I want a nXmXp matrix where n is condtion (of which we have 4),
% m is opponent (of which we have 3), and p is trial (of which we have a
% variable number, but I will cap it at 40 because that should be the
% average sincel that's the number of total trails divided by the numner of
% opponents)
all_subs = unique(soc(:, 7));
nsub = length(all_subs);

% subject
% social nonsocial
% opponent
% trial
ugm = NaN(nsub, 2, 3, ntr);
ugs = NaN(nsub, 2, 3, ntr);

% ugs = NaN(4, 3, ntr, length(unique(nonsoc(:, 7))));
for ss = 1:nsub
    for oo = 1:size(ugm, 3)
        %     for tt = 1:size(ugm, 3)
        
        tmp1 = (soc(soc(:, 5) == oo & soc(:, 7) == all_subs(ss), 1));
        tmp2 = (nonsoc(nonsoc(:, 5) == oo & nonsoc(:, 7) == all_subs(ss), 1));

        ugm(ss, 1, oo, 1:ntr) = tmp1(1:ntr);
        ugm(ss, 2, oo, 1:ntr) = tmp2(1:ntr);
        
        
       
%         ugm(ss, 2, oo, 1:ntr) = mean(nonsoc(nonsoc(:, 5) == oo & nonsoc(:, 7) == all_subs(ss) & nonsoc(:, 8) <= ntr, 1)); 
        
        %         ugm(3, oo, tt) = mean(soc_exp(soc_exp(:, 5) == oo & soc_exp(:, 8) == tt, 1));          % social explicit
        %         ugm(4, oo, tt) = mean(nonsoc_exp(nonsoc_exp(:, 5) == oo & nonsoc_exp(:, 8) == tt, 1)); % non social explicit
        %
        %         % standar errors
        %         ugs(1, oo, tt) = ste(soc_imp(soc_imp(:, 5) == oo & soc_imp(:, 8) == tt, 1));          % social implicit
        %         ugs(2, oo, tt) = ste(nonsoc_imp(nonsoc_imp(:, 5) == oo & nonsoc_imp(:, 8) == tt, 1)); % non social implicit
        %
        %         ugs(3, oo, tt) = ste(soc_exp(soc_exp(:, 5) == oo & soc_exp(:, 8) == tt, 1));          % social explicit
        %         ugs(4, oo, tt) = ste(nonsoc_exp(nonsoc_exp(:, 5) == oo & nonsoc_exp(:, 8) == tt, 1)); % non social explicit
        %
        
    end
    %     end
end

%% plots that work!
close all
all_leg = {'r-^', 'b-.o', 'k:s'};

figure
subplot(1, 2, 1)


hold on
for ii = 1:3
    
%     plot(1:ntr, mean(squeeze((ugm(:, 2, ii, :)))), all_leg{ii}, 'markers',5, 'linewidth', 1)
    
    errorbar(mean(squeeze((ugm(:, 2, ii, :)))), ste(squeeze((ugm(:, 2, ii, :)))), all_leg{ii}, 'markers',8, 'linewidth', 1)
    
end
ylim([4 10])
title('Non social')
hold off

subplot(1, 2, 2)
hold on
for ii = 1:3
  
%     plot(1:ntr, mean(squeeze((ugm(:, 1, ii, :)))), ste(squeeze((ugm(:, 1, ii, :)))), all_leg{ii}, 'markers',5, 'linewidth', 1)
    
%     errorbar(mean(squeeze((ugm(:, 1, ii, :)))), ste(squeeze((ugm(:, 1, ii, :)))), all_leg{ii}, 'markers',8, 'linewidth', 1)
    errorbar(mean(squeeze((ugm(:, 1, ii, :)))), ste(squeeze((ugm(:, 1, ii, :)))), all_leg{ii}, 'markers',8, 'linewidth', 1)
    
end
ylim([4 10])
title('Social')
hold off
%%
figure
for pp = 1:2
    subplot(2, 2, pp)
    hold on
    for ii = 1:3
%         errorbar(squeeze((ugm(pp, ii, :))), squeeze((ugs(pp, ii, :))), all_leg{ii}, 'markers',5, 'linewidth', 1)
        plot(1:my_x_lim, squeeze((snugs(pp, ii, :))), all_leg{ii}, 'markers',5, 'linewidth', 1)
    end
    legend(plot_leg)
    ylim([3 10])
    xlim([0 my_x_lim])
    hold off
    title(plot_titles{pp})
end
    
%% NICE PLOT OF ALL DATA
close all

my_x_lim = 40;


plot_titles = {'SocImp', 'NonSocImp', 'SocExp', 'NonSocExp'};
plot_leg = {'0 endow', '10 endow', '20 endow'};
all_leg = {'r-^', 'b-.o', 'k:s'};

figure
for pp = 1:4
    subplot(2, 2, pp)
    hold on
    for ii = 1:3
%         errorbar(squeeze((ugm(pp, ii, :))), squeeze((ugs(pp, ii, :))), all_leg{ii}, 'markers',5, 'linewidth', 1)
        plot(1:my_x_lim, squeeze((ugm(pp, ii, :))), all_leg{ii}, 'markers',5, 'linewidth', 1)
    end
    legend(plot_leg)
    ylim([3 10])
    xlim([0 my_x_lim])
    hold off
    title(plot_titles{pp})
end
    

%---------------%---------------%---------------%---------------%---------------
%---------------%---------------%---------------%---------------%---------------
%---------------%---------------%---------------%---------------%---------------
%---------------%---------------%---------------%---------------%---------------
%---------------%---------------%---------------%---------------%---------------
%%
y = [mean(soc_imp(soc_imp(:, 5) == 1, 1)),  mean(soc_imp(soc_imp(:, 5) == 2, 1)), mean(soc_imp(soc_imp(:, 5) == 3, 1));...
    mean(nonsoc_imp(nonsoc_imp(:, 5) == 1, 1)),  mean(nonsoc_imp(nonsoc_imp(:, 5) == 2, 1)), mean(nonsoc_imp(nonsoc_imp(:, 5) == 3, 1));...
    mean(soc_exp(soc_exp(:, 5) == 1, 1)),  mean(soc_exp(soc_exp(:, 5) == 2, 1)), mean(soc_exp(soc_exp(:, 5) == 3, 1));...
    mean(nonsoc_exp(nonsoc_exp(:, 5) == 1, 1)),  mean(nonsoc_exp(nonsoc_exp(:, 5) == 2, 1)), mean(nonsoc_exp(nonsoc_exp(:, 5) == 3, 1))];

stey = [ste(soc_imp(soc_imp(:, 5) == 1, 1)),  ste(soc_imp(soc_imp(:, 5) == 2, 1)), ste(soc_imp(soc_imp(:, 5) == 3, 1));...
    ste(nonsoc_imp(nonsoc_imp(:, 5) == 1, 1)),  ste(nonsoc_imp(nonsoc_imp(:, 5) == 2, 1)), ste(nonsoc_imp(nonsoc_imp(:, 5) == 3, 1));...
    ste(soc_exp(soc_exp(:, 5) == 1, 1)),  ste(soc_exp(soc_exp(:, 5) == 2, 1)), ste(soc_exp(soc_exp(:, 5) == 3, 1));...
    ste(nonsoc_exp(nonsoc_exp(:, 5) == 1, 1)),  ste(nonsoc_exp(nonsoc_exp(:, 5) == 2, 1)), ste(nonsoc_exp(nonsoc_exp(:, 5) == 3, 1))];


% bar(y)
% hold on
leg_txt = {'soc imp', 'nonsoc imp', 'soc exp', 'nonsoc exp'};
point_spacing = [.8 1.8 2.8];
all_leg = {'r-s', 'b-o', 'r:s', 'b:o'};
hold on
for ii = 1:4
    errorbar(point_spacing, y(ii, :), stey(ii, :), all_leg{ii}, 'markers',12, 'linewidth', 2)
    point_spacing = point_spacing + .1;
end
legend(leg_txt)
set(gca, 'XTick', [1:3], 'fontsize', 20)
hold off

%%






for k_sub = 1:nsub
    
    %% check sub
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        %% Get UG data
        flnm = fullfile(data_dir,strcat('Sub_',sub_nm,'_UG.mat'));
        load(flnm)
        
        
        
    end
end
        
        %%
        
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

Gsub = ~isnan(good_sub);
nGsub = sum(double(Gsub));

save('Learning_All_2017_12_08')