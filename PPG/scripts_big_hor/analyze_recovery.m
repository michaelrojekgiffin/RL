clear

clc


load 'recovery60sims_2018_12_03.mat'

close all force

% for predator
for k_true = modelspace
    
    MP = [Px_rnd,Pa_rnd,Pb_rnd,Plr1_rnd,Plr2_rnd];
    MP = MP(1:length(MP)/2, :); % only first half is for predators
    
    % laplace
    LL = squeeze(LPP(:,k_true,:, 1)); % 1st page of 4th dimension is predator
    
    % non-laplace
%     LL = squeeze(ll(:,k_true,:, 1)); % 1st page of 4th dimension is predator
    LL(any(isnan(LL), 2), :) = [];
    
    for k_est= modelspace
        bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
    end
    
%     [postBMC,outBMC]=VBA_groupBMC(-LL');
    [postBMC,outBMC]=VBA_groupBMC(-bic');
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
c.Label.String = 'Exceedance probability Predator (%)';




%% for prey
for k_true = modelspace
    
    MP = [Px_rnd,Pa_rnd,Pb_rnd,Plr1_rnd,Plr2_rnd];
    MP = MP((length(MP)/2)+1:end, :); % second half for prey
    
    % laplace
    LL = squeeze(LPP(:,k_true,:, 2)); % 2nd page of 4th dimension is predator
    
    % non-laplace
%     LL = squeeze(ll(:,k_true,:, 2)); % 2nd page of 4th dimension is predator
    LL(any(isnan(LL), 2), :) = [];
    
    for k_est= modelspace
        bic(:,k_est)=-2*-LL(:,k_est) + nfpm(k_est)*log(nc*n_trial*n_sess); % l2 is already positive
    end
    
%     [postBMC,outBMC]=VBA_groupBMC(-LL');
    [postBMC,outBMC]=VBA_groupBMC(-bic');
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
c.Label.String = 'Exceedance probability Prey (%)';


%% Predator parameters
close all
% param 

MP = [Px_rnd,Pa_rnd,Pb_rnd,Plr1_rnd,Plr2_rnd];
MP = MP(1:length(MP)/2, :); % only first half is for predators

for k_true = modelspace
    
    
    legB = {'rating temperature','Sig intercept','Sig Slope','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    
    
    title(strcat(['Model ',num2str(k_true)]));
    for k = 1:5
        tmp_params = squeeze(parameters(:,k,k_true,k_true, 1));
        tmp_params(any(isnan(tmp_params), 2), :) = [];
        
        tmp_LPP = squeeze(parametersLPP(:,k,k_true,k_true, 1));
        tmp_LPP(any(isnan(tmp_LPP), 2), :) = [];
        
        subplot(2,5,k)
        plot(MP(:,k),tmp_params,'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(MP(:,k),tmp_params);
        txt1 = sprintf('r = %f\n p = %f', corrR(k), corrP(k));
        text(mean(MP(:,k)), mean(tmp_params), txt1);
        lsline;
        if k == 1
            title(['Predator Model ' num2str(k_true)]);
        end
        
        
        subplot(2,5,5+k)
        plot(MP(:,k) ,tmp_LPP,'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),tmp_LPP);
        txt1 = sprintf('r = %f\n p = %f', corrR_LPP(k), corrP_LPP(k));
        text(mean(MP(:,k)), mean(tmp_LPP), txt1);
        lsline;
        
        
    end
    
end



%% Prey parameters

% param 

MP = [Px_rnd,Pa_rnd,Pb_rnd,Plr1_rnd,Plr2_rnd];
MP = MP((length(MP)/2)+1:end, :); % second half for prey

for k_true = modelspace
    
    
    legB = {'rating temperature','Sig intercept','Sig Slope','learning rate 1','learning rate 2'};
    
    figure;
    set(gcf,'Color',[1,1,1])
    
    
    title(strcat(['Model ',num2str(k_true)]));
    for k = 1:5
        tmp_params = squeeze(parameters(:,k,k_true,k_true, 2)); % last dimension is 2 for prey
        tmp_params(any(isnan(tmp_params), 2), :) = [];
        
        tmp_LPP = squeeze(parametersLPP(:,k,k_true,k_true, 2)); % last dimension is 2 for prey
        tmp_LPP(any(isnan(tmp_LPP), 2), :) = [];
        
        subplot(2,5,k)
        plot(MP(:,k),tmp_params,'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k}]));
        [corrR(k),corrP(k)] = corr(MP(:,k),tmp_params);
        txt1 = sprintf('r = %f\n p = %f', corrR(k), corrP(k));
        text(mean(MP(:,k)), mean(tmp_params), txt1);
        lsline;
        if k == 1
            title(['Prey Model ' num2str(k_true)]);
        end
        
        
        subplot(2,5,5+k)
        plot(MP(:,k) ,tmp_LPP,'o',...
            'MarkerEdgeColor',[0,0,0],...
            'MarkerFaceColor',[1,1,1])
        xlabel(strcat(['true ' legB{k}]));
        ylabel(strcat(['estimated ' legB{k} ' LPP']));
        
        [corrR_LPP(k),corrP_LPP(k)] = corr(MP(:,k),tmp_LPP);
        txt1 = sprintf('r = %f\n p = %f', corrR_LPP(k), corrP_LPP(k));
        text(mean(MP(:,k)), mean(tmp_LPP), txt1);
        lsline;
        
        
    end
    
end