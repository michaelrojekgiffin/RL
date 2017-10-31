clc
clear
close all

rng('shuffle')
options         = optimset('Algorithm', 'interior-point', 'MaxIter', 1000000);
% Params
nsub   = 40;
% PARAMS = [10*rand(nsub,1)-10,2*rand(nsub,1),1+3*rand(nsub,1)];
PARAMS = [10*rand(nsub,1)-10,2*rand(nsub,1),1+3*rand(nsub,1)];

tr_choice = repmat(nchoosek(0:1:20,2),1,1);


% model specs
offers  = 0:1:20;
endow   = 20*ones(1,numel(offers));
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% pre-allocate
% sub_o_freq  = zeros(nsub, 11);
kO          = NaN(nsub,numel(offers));
PA_sub      = NaN(nsub,numel(offers));
EV_sub      = NaN(nsub,numel(offers));
sub_o       = NaN(nsub,size(tr_choice,1));
pc          = NaN(nsub,numel(offers));

parameters      = NaN(nsub,3);  parametersLPP   = NaN(nsub,3);
ll              = NaN(nsub,1);  LPP             = NaN(nsub,1);


%subject loop
for k_sub   = 1:nsub
    
    % get the distribs
    PA     = logitp([PARAMS(k_sub,1),PARAMS(k_sub,2)],offers);            % compute proba of accepting the offers given current model
    EV     = (endow - offers).* PA ;                                   % compute EV of the offers given current model
    
    PA_sub(k_sub,:) = PA;
    EV_sub(k_sub,:) = EV;
    
    for k_tr= 1:size(tr_choice,1)
        pc(k_sub,k_tr)     = exp(PARAMS(k_sub,3).*EV(tr_choice(k_tr,1)+1) ) ./ sum(exp(PARAMS(k_sub,3).*EV(tr_choice(k_tr,:)+1)));
        sub_o(k_sub,k_tr)   = rand()<pc(k_sub,k_tr) ;                                   % resample Offer in pdf (="soft-max")
    end
    
    % fit the prior
    
    n_rep           = 5;
    parameters_rep  = NaN(n_rep,3);     parametersLPP_rep  = NaN(n_rep,3);
    ll_rep          = NaN(n_rep,1);     LPP_rep             = NaN(n_rep,1);
    
    lb = [-20 0 0];     LB = [-Inf 0 0];
    ub = [10 20 20];    UB = [Inf Inf Inf];
    
        
    ub = [10 20 20];
    
    for k_rep = 1:n_rep
%         x0 = [-5*rand() 5*rand()  5*rand()];
        x0 = [-5*rand() rand()  5*rand()];
        
        % standard estim
        [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) estimate_priors_bin(x,tr_choice,sub_o(k_sub,:)),x0,[],[],[],[],lb,ub,[],options);
        % laplace approximation
        [parametersLPP_rep(k_rep,1:3),LPP_rep(k_rep,1)]=fmincon(@(x) laplace_priors(x,tr_choice,sub_o(k_sub,:)),x0,[],[],[],[],LB,UB,[],options);
    end
    [~,pos]                  = min(ll_rep);
    parameters(k_sub,:)      = parameters_rep(pos(1),:);
    ll(k_sub)                = ll_rep(pos(1),:);
    
    [~,posLPP]                  = min(LPP_rep);
    parametersLPP(k_sub,:)      = parametersLPP_rep(posLPP(1),:);
    LPP(k_sub)                 = LPP_rep(posLPP(1),:);
    
end


figure;
set(gcf,'Color',[1,1,1])
% subplot(1,3,1)
% hold on
% mtp = mean(sub_o_freq);
% stp = std(sub_o_freq)./sqrt(nsub);
% bar(offers,mtp,'FaceColor',.7.*[1,1,1])
% errorbar(offers,mtp,stp,'k','LineStyle','none')
% plot(offers,mean(pc),'-o',...
%     'Color',0.*[1,1,1],...
%     'MarkerFaceColor',1.*[1,1,1],...
%     'MarkerEdgeColor',0.*[1,1,1])
% set(gca,'XLim',[-.5 10.5],...
%     'XTick',0:10)
% xlabel('offers')
% ylabel('frequency (%)')
% legend('Emprirical','Model')
%
subplot(1,3,2)
plot(offers,mean(EV_sub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated Expected value')

subplot(1,3,3)
plot(offers,mean(PA_sub),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')



figure;
param_names = {'intercept', 'slope', 'error'};
set(gcf,'Color',[1,1,1])
for k = 1:3
    subplot(2,3,k)
    plot(PARAMS(:,k),parameters(:,k),'o',...
        'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0,0,0])
    [corrR(k),corrP(k)] = corr(PARAMS(:,k),parameters(:,k));
    txt1 = sprintf('r = %f\n p = %f', corrR(k),corrP(k));
    title(txt1);
    lsline;
    xlabel(strcat(['true ' param_names{k}]));
    ylabel(strcat(['estimated ' param_names{k}]));
    
        
    subplot(2,3,3+k)
    plot(PARAMS(:,k),parametersLPP(:,k),'o',...
        'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0,0,0])
    [corrR_lpp(k),corrP_lpp(k)] = corr(PARAMS(:,k),parametersLPP(:,k));
    txt1 = sprintf('r = %f\n p = %f', corrR_lpp(k),corrP_lpp(k));
    title(txt1);
    lsline;
    xlabel(strcat(['true ' param_names{k}]));
    ylabel(strcat(['estimated ' param_names{k} ' LPP']));
    
    
end


