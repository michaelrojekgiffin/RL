% to be run from the directory where this script is stored
% this script reads in the subject's data (priors and UG task), fits our 4
% models to the data, and outputs some plots and stuff
clear
close all
clc
% randomize generator seed
%--------------------------
rng('shuffle')
options         = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000, 'display', 'off');
%------------------------------------------
% change depending on where the subject's data is stored
%------------------------------------------
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '3cond_interleafed' filesep 'processed' ];
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '3cond_blocked' filesep 'processed' ];
data_path = ['..' filesep 'data' filesep 'processed' ];
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '60_trials' filesep '3cond_blocked' filesep 'processed' ];
%------------------------------------------
%------------------------------------------

priors_soc_dir    = dir([data_path filesep '*_SocPriors*']);
priors_nonsoc_dir = dir([data_path filesep '*NonSocPriors*']);
UG_dir            = dir([data_path filesep '*UG*']);

modelspace = [1 2 3 4];
nfpm=[2 3 2 3];
nmods = numel(modelspace);

for ss = 1:length(UG_dir)
    sub_name = UG_dir(ss).name(1:12);
    
    all_there = 0;
    for fsub = 1:length(priors_soc_dir)
        tmpidx = (strcmp(priors_soc_dir(fsub).name(1:12), sub_name));
        if tmpidx == 1
            soc_priors_idx = tmpidx;
            all_there = all_there+1;
            
            priors_soc_IDX = fsub;
        end
    end
    for fsub = 1:length(priors_nonsoc_dir)
        tmpidx = (strcmp(priors_nonsoc_dir(fsub).name(1:12), sub_name));
        if tmpidx == 1
            nonsoc_priors_idx = tmpidx;
            all_there = all_there+1;
            
            priors_nonsoc_IDX = fsub;
        end
    end
    % only run the following if the subject has both priors and UG data
    if all_there == 2
        
        % these are the data files, for the priors these are what's contained
        % column 1: option 1
        % column 2: option 2
        % column 3: offer
        % column 4: accepted (1 for yes, 0 for no)
        % column 5: payment
        load(fullfile(data_path, priors_soc_dir(priors_soc_IDX).name));
        priors_soc_data = subdata;
        if ~isnan(sub_age(1))
            rl_sub_age = sub_age(1);
        end
        
        load(fullfile(data_path, priors_nonsoc_dir(priors_nonsoc_IDX).name));
        priors_nonsoc_data = subdata;
        if ~isnan(sub_age(1))
            rl_sub_age = sub_age(1);
        end
        
        % only run it if the subject is within our age limits, 18-25
        if rl_sub_age < 34
            
            % for the UG datafile
            % column 1: offer
            % column 2: accepted (1 for yes, 0 for no)
            % column 3: payment
            % column 4: social (1 for social, 0 for non-social)
            % column 5: opponent, 1 for starting endowment of 0, 2 for starting
            % endowment of 10, and 3 for starting endowment of 20
            % column 6: explicit (1 for explicit, 0 for implicit)
            load(fullfile(data_path, UG_dir(ss).name));
            UG_data = subdata;
            
            
            %--------------------------------------------------------------- 
            if UG_data(1, 6) == 0
            
                
                
            
            
            sub_name = UG_dir(ss).name(1:12);
            
            soc   = UG_data(UG_data(:, 4) == 1, :);
            soc_o = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            soc_r = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            soc_p = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            for soso = 1:size(soc_o, 2)
                soc_o(1:length(soc(soc(:, 5) == soso, 1)), soso) = soc(soc(:, 5) == soso, 1);
                soc_r(1:length(soc(soc(:, 5) == soso, 1)), soso) = soc(soc(:, 5) == soso, 2);
                soc_p(1:length(soc(soc(:, 5) == soso, 1)), soso) = soc(soc(:, 5) == soso, 3);
            end
            all_soc_o(:, :, ss) = soc_o;
            
            
            nonsoc   = UG_data(UG_data(:, 4) == 0, :);
            nonsoc_o = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            nonsoc_r = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            nonsoc_p = NaN(length(UG_data), length(unique(UG_data(:, 5))));
            for soso = 1:size(soc_o, 2)
                nonsoc_o(1:length(nonsoc(nonsoc(:, 5) == soso, 1)), soso) = nonsoc(nonsoc(:, 5) == soso, 1);
                nonsoc_r(1:length(nonsoc(nonsoc(:, 5) == soso, 1)), soso) = nonsoc(nonsoc(:, 5) == soso, 2);
                nonsoc_p(1:length(nonsoc(nonsoc(:, 5) == soso, 1)), soso) = nonsoc(nonsoc(:, 5) == soso, 3);
            end
            all_nonsoc_o(:, :, ss) = nonsoc_o;
            
            soc_priors_o = NaN(1, length(priors_soc_data));
            for ii = 1:length(priors_soc_data)
                if priors_soc_data(ii,3) == priors_soc_data(ii,1)
                    soc_priors_o(ii) = 0;
                else
                    soc_priors_o(ii) = 1;
                end
            end
            nonsoc_priors_o = NaN(1, length(priors_nonsoc_data));
            for ii = 1:length(priors_nonsoc_data)
                if priors_nonsoc_data(ii,3) == priors_nonsoc_data(ii,1)
                    nonsoc_priors_o(ii) = 0;
                else
                    nonsoc_priors_o(ii) = 1;
                end
            end
            
            % fit the prior
            n_rep                    = 5;
            parametersLPP_rep_soc    = NaN(n_rep,3);
            parametersLPP_rep_nonsoc = NaN(n_rep,3);
            LPP_rep_soc              = NaN(n_rep,1);
            LPP_rep_nonsoc           = NaN(n_rep,1);
            
            LB = [-Inf 0 0];
            UB = [Inf Inf Inf];
            
            
            for k_rep = 1:n_rep
                x0 = [-5*rand() rand()  5*rand()];
                
                % laplace approximation
                [parametersLPP_rep_soc(k_rep,1:3),LPP_rep_soc(k_rep,1)]=fmincon(@(x) laplace_priors_priors(x,priors_soc_data(:,1:2),soc_priors_o),x0,[],[],[],[],LB,UB,[],options);
                
                [parametersLPP_rep_nonsoc(k_rep,1:3),LPP_rep_nonsoc(k_rep,1)]=fmincon(@(x) laplace_priors_priors(x,priors_nonsoc_data(:,1:2),nonsoc_priors_o),x0,[],[],[],[],LB,UB,[],options);
            end
            
            [~,posLPP]            = min(LPP_rep_soc);
            soc_priors(ss,:)      = parametersLPP_rep_soc(posLPP(1),1:2);
            soc_priors_LPP(ss)    = LPP_rep_soc(posLPP(1),:);
            
            [~,posLPP]            = min(LPP_rep_nonsoc);
            nonsoc_priors(ss,:)   = parametersLPP_rep_nonsoc(posLPP(1),1:2);
            nonsoc_priors_LPP(ss) = LPP_rep_nonsoc(posLPP(1),:);
            
            
            for kem = modelspace
                fprintf('fitting model %d for subject %s, %d out of %d\n', kem, sub_name, ss, length(UG_dir));
                n_rep              = 5;
                soc_parametersLPP_rep  = NaN(n_rep,3);
                socLPP_rep            = NaN(n_rep,1);
                nonsoc_parametersLPP_rep  = NaN(n_rep,3);
                nonsocLPP_rep            = NaN(n_rep,1);
                
                LB = [0 0 0];
                UB = [Inf 1 1];
                lb = [0 0 0];
                ub = [15 1 1];
                ddb = ub - lb;
                
                
                for k_rep = 1:n_rep
                    x0 = lb + rand(1,3).*ddb;
                    % %laplace estimation
                    % the function takes as input:
                    % parameters
                    % the offers made by the subjects
                    % whether it was accepted or not
                    % prior intercept
                    % prior slope
                    % model number
                    [soc_parametersLPP_rep(k_rep,1:3),socLPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2(x,soc_o,soc_r,soc_priors(ss,1),soc_priors(ss,2),kem),x0,[],[],[],[],LB,UB,[],options);
                    
                    [nonsoc_parametersLPP_rep(k_rep,1:3),nonsocLPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors_learning2(x,nonsoc_o,nonsoc_r,nonsoc_priors(ss,1),nonsoc_priors(ss,2),kem),x0,[],[],[],[],LB,UB,[],options);
                end
                
                [~,socposLPP]                       = min(socLPP_rep);
                soc_parametersLPP(ss,:,kem)      =   soc_parametersLPP_rep(socposLPP(1),:);
                socLPP(ss,kem)                      =   socLPP_rep(socposLPP(1),:);
                
                [~,nonsocposLPP]                       = min(nonsocLPP_rep);
                nonsoc_parametersLPP(ss,:,kem)      =   nonsoc_parametersLPP_rep(nonsocposLPP(1),:);
                nonsocLPP(ss,kem)                      =   nonsocLPP_rep(nonsocposLPP(1),:);
            end
            
            
            %--------------------------------------------------------------- 
            end 
        end
    end
end

%%
options.families = {[1,2], [3,4]} ;

rowsWithZeros = any(socLPP==0, 2);
socLPP     = socLPP(~rowsWithZeros, :);

LL = socLPP;

for k_est= modelspace
    [~, bic(:,k_est)] = aicbic(LL(:, k_est), nfpm(k_est), 120);
end

[postBMC,outBMC]=VBA_groupBMC(-bic'./2);


rowsWithZeros = any(nonsocLPP==0, 2);
nonsocLPP     = nonsocLPP(~rowsWithZeros, :);

LL = nonsocLPP;

for k_est= modelspace
    [~, bic(:,k_est)] = aicbic(LL(:, k_est), nfpm(k_est), 120);
end

[postBMC,outBMC]=VBA_groupBMC(-bic'./2);


% % [postBMC,outBMC]=VBA_groupBMC(-LL');
% BMC_output(k_true).post = postBMC;
% BMC_output(k_true).out = outBMC;
%
% Ep(k_true,:) = 100*BMC_output(k_true).out.ep;

%%
% close all
% 
rowsWithZeros = any(soc_parametersLPP(:,2,3)==0, 2);
soc_parametersLPP3     = soc_parametersLPP(~rowsWithZeros, 2, 3);

rowsWithZeros = any(nonsoc_parametersLPP(:,2,3)==0, 2);
nonsoc_parametersLPP3     = nonsoc_parametersLPP(~rowsWithZeros, 2, 3);

mean(soc_parametersLPP3)
mean(nonsoc_parametersLPP3)
% dlmwrite('~/Desktop/learningrates.txt', [soc_parametersLPP(:,2,3), nonsoc_parametersLPP(:,2,3)], 'delimiter','\t');

[H, P ci stats] = ttest(soc_parametersLPP3, nonsoc_parametersLPP3, 'Tail','left')
%%

% [H, P] = ttest(soc_parametersLPP(:,1,3), nonsoc_parametersLPP(:,1,3))
mlr3 = [mean(soc_parametersLPP(:,2,3)), mean(nonsoc_parametersLPP(:,2,3))];
selr3 = [(std(soc_parametersLPP(:,2,3)))/(sqrt(length(soc_parametersLPP(:,2,3)))), (std(nonsoc_parametersLPP(:,2,3)))/(sqrt(length(nonsoc_parametersLPP(:,2,3))))];

figure
hold on;
bar(1:2, mlr3)
errorbar(1:2, mlr3, selr3, '.k', 'linewidth', 2)

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'Social' 'Nonsocial'})
txt = sprintf('p = %.4f', P);
text(1.2, .12, txt, 'Color', 'r', 'fontsize', 12);

xlabel('Condition')
ylabel('Learning rate')
title('Learning rate in model 3')
hold off;

%%
cond1 = (squeeze(all_nonsoc_o(:, 1, :)));
cond2 = (squeeze(all_nonsoc_o(:, 2, :)));
cond3 = (squeeze(all_nonsoc_o(:, 3, :)));

for ii = 1:size(cond1, 2)
    tc1 = cond1;
    tc1(any(isnan(tc1), 2), :) = [];
    mc1(ii) = mean(tc1(ii));
    
    tc2 = cond2;
    tc2(any(isnan(tc2), 2), :) = [];
    mc2(ii) = mean(tc2(ii));
    
    tc3 = cond3;
    tc3(any(isnan(tc3), 2), :) = [];
    mc3(ii) = mean(tc3(ii));
end

nmc1 = mean(mc1);
nmc2 = mean(mc2);
nmc3 = mean(mc3);


cond1 = (squeeze(all_soc_o(:, 1, :)));
cond2 = (squeeze(all_soc_o(:, 2, :)));
cond3 = (squeeze(all_soc_o(:, 3, :)));

for ii = 1:size(cond1, 2)
    tc1 = cond1;
    tc1(any(isnan(tc1), 2), :) = [];
    mc1(ii) = mean(tc1(ii));
    
    tc2 = cond2;
    tc2(any(isnan(tc2), 2), :) = [];
    mc2(ii) = mean(tc2(ii));
    
    tc3 = cond3;
    tc3(any(isnan(tc3), 2), :) = [];
    mc3(ii) = mean(tc3(ii));
end

smc1 = mean(mc1);
smc2 = mean(mc2);
smc3 = mean(mc3);


fprintf('social cond 1: %d\nsocial cond 2: %f\nsocial cond 3: %f\nnonsocial cond1: %f\nnonsocial cond 2: %f\nnonsocial cond 3: %f\n', smc1, smc2, smc3, nmc1, nmc2, nmc3);