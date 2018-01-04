clear

clc


load 'ML_recovery_2017_10_19.mat'

close all force

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


figure
set(gcf,'Color',[1,1,1])


colormap(flipud(gray))
imagesc(flipud(Ep))
ylabel('simulated model #')
xlabel('estimated model #')
set(gca,'XTick',1:4,...
    'YTick',1:4,...
    'XTickLabel',(1:4),...
    'YTickLabel',fliplr(1:4))

c = colorbar;
c.Label.String = 'Exceedance probability (%)';




% param 

for k_true = modelspace
    
    
    legB = {'rating temperature','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    
    
    title(strcat(['Model ',num2str(k_true)]));
    for k = 1:3
        
        subplot(2,3,k)
        plot(MP(:,k),squeeze(parameters(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(MP(:,k),squeeze(parameters(:,k,k_true,k_true)));
        
        subplot(2,3,3+k)
        plot(MP(:,k) ,squeeze(parametersLPP(:,k,k_true,k_true)),'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),squeeze(parametersLPP(:,k,k_true,k_true)));
        
    end
    
end
