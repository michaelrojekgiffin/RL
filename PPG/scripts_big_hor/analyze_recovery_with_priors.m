close all
clear all
clc

% analyzes the output of test_recovery_with_priors.m
project_name = '~/Dropbox/RL/PPG/'; % for use in michael's dropbox


% load([project_name 'arxiv/test_recovery_workspace_07-Nov-2017.mat'])
load(['test_recovery_workspace_30_14-Mar-2018.mat'])

%%

% Recovery plot Loop
%----------
for predprey_count = 1:length(predprey_array)
    predprey = predprey_array{predprey_count};
  
    
    for model_iter = 1:length(nmodel_array)
        nmodel = nmodel_array(model_iter);
        
        % since the different models have different number of parameters we
        % have to specify this here
        if nmodel == 1 || nmodel == 3
            numfreeparams = 2;
            PP = [Px_rnd,Plr1_rnd];         % true parameters
            legB = {'rating temperature','learning rate 1'};
        elseif nmodel == 2 || nmodel == 4
            numfreeparams = 3;
            PP = [Px_rnd,Plr1_rnd,Plr2_rnd];         % true parameters
            legB = {'rating temperature','learning rate 1','learning rate 2'};
        end
        
        
        rec_plot = figure;
        set(gcf,'Color',[1,1,1])

         % Make the subplots
        for k = 1:numfreeparams
            
            subplot(2,numfreeparams,k)
            plot(PP(:,k),con_parameters(predprey_count, :, k, nmodel, nmodel),'o',... % correlate across all simulations for specific parameter 
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,1])
            xlabel(strcat(['true ' legB{k}]));
            ylabel(strcat(['estimated ' legB{k}]));
            [corrR(k),corrP(k)] = corr(PP(:,k),con_parameters(predprey_count, :, k, nmodel, nmodel)');
            
            txt1 = sprintf('r = %f\n p = %f', corrR(k),corrP(k));
            if k == 1
                title(sprintf('%s model %d\n %s', predprey, nmodel, txt1));
            else
                title(txt1)
            end
            %     text(mean(PP(:,k)), mean(parametersLPP(:,k)), txt1);
            lsline;
            
            subplot(2,numfreeparams,numfreeparams+k)
            plot(PP(:,k) ,con_parametersLPP(predprey_count, :, k, nmodel, nmodel),'o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,1])
            xlabel(strcat(['true ' legB{k}]));
            ylabel(strcat(['estimated ' legB{k} ' LPP']));
            
            [corrR_LPP(k),corrP_LPP(k)] = corr(PP(:,k),con_parametersLPP(predprey_count, :, k, nmodel, nmodel)');
            txt1 = sprintf('r = %f\n p = %f', corrR_LPP(k),corrP_LPP(k));
            title(txt1)
            %     text(mean(PP(:,k)), mean(parametersLPP(:,k)), txt1);
            lsline;
            
        end
        
            
    end
end

%% Confusion matrix
% confusion matrix
% parameters and parameresLPP are each mXnXpXz matrices have the following
% dimensions 
% dim 1: simulated subject (so each entry is one simulated subject)
% dim 2: paraemters - so each entry is a different parameter of the model
% dim 3: the actual model we used to generate the simulated data 
% dim 4: the parameters estimated from 
%----------------%----------------%----------------%----------------%----------------
% thereofer, parameters(1, 1, 1, 1) will give me parameter 1 estimated from
% the simulated data assuming the data were generated using model 1, when
% the data were actually generated from model 1, for parameter 1 (beta), of
% simulated subject 1. parameters(1, 1, 1, 2) will give me a parameter 1
% assuming the data were generated using model 2 when data were really
% generated using model 1 of simulated subject 1. 
%
% A similar story is true of the variable ll. ll(1, 1, 1) gives me the
% log-liklihood of the data for simulated subject 1 when their data were
% generated using model 1 and we are assuming their data were generated
% using model 1. ll(1, 1, 2) gives us the log-liklihood of subject 1 for
% data generated using model 1 but assuming it was generated using model 2.
%----------------%----------------%----------------%----------------%----------------


BIC = NaN(n_sims, 4, 4);
for ppg = 1:length(predprey_array)
    predprey = predprey_array{ppg};
    for k_true = 1:4
        nfpm=[2 3 2 3]; % number of free parameters
        for k_est=1:4
            [~, BIC(:,k_true, k_est)] = aicbic(-con_LPP(ppg, :, k_true, k_est), nfpm(k_est), n_trial);
            
% %             bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
            
            BIC_mean(k_true, k_est)          = mean(BIC(:, k_true, k_est));
            ll_mean(k_true, k_est)           = mean(-con_LPP(ppg, :, k_true, k_est));
            
            BIC_sum(k_true, k_est)          = sum(BIC(:, k_true, k_est));

        end
    end
    figure
    imagesc(BIC_sum)
    colsize = size(BIC_sum, 2);
    colcounter = 0;
    rowcounter = 1;
    for ii = 1:numel(BIC_sum)
        colcounter = colcounter+1;
        text(colcounter, rowcounter, num2str(BIC_sum(rowcounter, colcounter)), 'Color', 'r')
        if colcounter == colsize
            colcounter = 0;
            rowcounter = rowcounter + 1;
        end
    end
    colormap((gray))
    title(sprintf('BIC %s', predprey))
    colorbar
    
    nfpm=[2 3 2 3]; % number of free parameters
    bic = NaN(n_sims, length(nmodel_array));
    for k_true = 1:length(nmodel_array)
        
        %     MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
        LL = squeeze(con_LPP(ppg, :,k_true,:));
%         LL = squeeze(con_ll(ppg, :,k_true,:));
        
        for k_est= 1:length(nmodel_array)
            bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
        end
        [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
        % [postBMC,outBMC]=VBA_groupBMC(-LL');
        BMC_output(k_true).post = postBMC;
        BMC_output(k_true).out = outBMC;
        
        Ep(k_true,:) = 100*BMC_output(k_true).out.ep;
    end
    
    figure
    set(gcf,'Color',[1,1,1])
    
    
    colormap(flipud(gray))
    imagesc(flipud(Ep))
    title(sprintf('confusion matrix %s', predprey))
    ylabel('simulated model #')
    xlabel('estimated model #')
    set(gca,'XTick',1:4,...
        'YTick',1:4,...
        'XTickLabel',(1:4),...
        'YTickLabel',fliplr(1:4))
    
    c = colorbar;
    c.Label.String = 'Exceedance probability (%)';
    

end

% ws_nm = ['test_recovery_workspace_' date];
% save(ws_nm)
