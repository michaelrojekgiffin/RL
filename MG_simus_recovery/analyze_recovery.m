clear
close all
clc

cur_dir = pwd;
all_data = dir([cur_dir, '/*MG_recovery*']);

for aa = 1:3
    
    flnm     = fullfile(cur_dir,all_data(aa).name);
    load(flnm)
    
    
    for k_true = modelspace
        
        MP = [Px_rnd,Plr1_rnd,Plr2_rnd];
        LL = squeeze(ll(:,k_true,:));
        
        for k_est= modelspace
            bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
        end
        [postBMC,outBMC]=VBA_groupBMC(-bic'./2);
        % [postBMC,outBMC]=VBA_groupBMC(-LL');
        BMC_output(k_true).post = postBMC;
        BMC_output(k_true).out = outBMC;
        
        Ep(k_true,:) = 100*BMC_output(k_true).out.ep;
        
    end
    
    
    conf_mat = figure;
    set(gcf,'Color',[1,1,1])
    
    
    colormap(flipud(gray))
    imagesc(flipud(Ep))
    ylabel('simulated model #')
    xlabel('estimated model #')
    set(gca,'XTick',1:4,...
        'YTick',1:4,...
        'XTickLabel',(1:4),...
        'YTickLabel',fliplr(1:4))
    title(sprintf('%d conditions', aa+2)); % +2 because conditions go from 3 to 5
    
    c = colorbar;
    c.Label.String = 'Exceedance probability (%)';
    
    print(conf_mat, ['..' filesep 'reports' filesep 'figures' filesep 'pres_2017_10_26' filesep 'conf_', num2str(aa+2), '_cond'], '-dpng');
    
    % param
    
    for k_true = modelspace
        
        
        legB = {'rating temperature','learning rate 1','learning rate 2'};
        
        rec_plot = figure;
        set(gcf,'Color',[1,1,1])
        
        
        title(strcat(['Model ',num2str(k_true)]));
        
        if k_true == 1 || k_true == 3
            numfreeparams = 2;
        else
            numfreeparams = 3;
        end
        
        for k = 1:numfreeparams
            
            subplot(2,numfreeparams,k)
            plot(MP(:,k),squeeze(parameters(:,k,k_true,k_true)),'o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,1])
            xlabel(strcat(['true ' legB{k}]));
            ylabel(strcat(['estimated ' legB{k}]));
            [corrR(k),corrP(k)] = corr(MP(:,k),squeeze(parameters(:,k,k_true,k_true)));
            
            txt1 = sprintf('r = %f\n p = %f', corrR(k),corrP(k));
            if k == 1
                title(sprintf('model %d\n%d environments\n %s', k_true, aa+2, txt1));
            else
                title(txt1)
            end
            lsline;
            
            subplot(2,numfreeparams,numfreeparams+k)
            plot(MP(:,k) ,squeeze(parametersLPP(:,k,k_true,k_true)),'o',...
                'MarkerEdgeColor',[0,0,0],...
                'MarkerFaceColor',[1,1,1])
            xlabel(strcat(['true ' legB{k}]));
            ylabel(strcat(['estimated ' legB{k} ' LPP']));
            
            [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),squeeze(parametersLPP(:,k,k_true,k_true)));
            
            txt1 = sprintf('r = %f\n p = %f', corrR_LPP(k),corrP_LPP(k));
            
            title(txt1)
            
            lsline;
            
        end
        print(rec_plot, ['..' filesep 'reports' filesep 'figures' filesep 'pres_2017_10_26' filesep 'recovery_model', num2str(k_true), '_',  num2str(aa+2), 'conds'], '-dpng');
        
    end
    
end