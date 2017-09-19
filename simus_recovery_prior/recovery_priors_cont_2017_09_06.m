clc
clear
close all

rng('shuffle')
options         = optimset('Algorithm', 'interior-point', 'MaxIter', 10000);
    
% Params
nsub   = 40;
ntr    = 10000;
PARAMS = [10*rand(nsub,1)-15,2*rand(nsub,1),1+3*rand(nsub,1)];

% model specs
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% pre-allocate
sub_o_freq  = zeros(nsub, 11);
kO          = NaN(nsub,ntr);
sub_o       = NaN(nsub,ntr);
PA          = NaN(nsub,numel(offers));
EV          = NaN(nsub,numel(offers));
pc          = NaN(nsub,numel(offers));

parameters  = NaN(nsub,3);
ll          = NaN(nsub,1);


%subject loop
for k_sub   = 1:nsub
    
    % get the distribs
    PA(k_sub,:)     = logitp([PARAMS(k_sub,1),PARAMS(k_sub,2)],offers);            % compute proba of accepting the offers given current model
    EV(k_sub,:)     = (endow - offers).* PA(k_sub,:) ;                                   % compute EV of the offers given current model
    pc(k_sub,:)     = exp(PARAMS(k_sub,3).*EV(k_sub,:) ) ./ sum(exp(PARAMS(k_sub,3).*EV(k_sub,:) ));
    
    pd              = makedist('multinomial','probabilities',pc(k_sub,:));         % estimate the pdf from the pC
    ypdf            = pdf(pd,1:numel(offers));                         % generate pdf for the offers
    kO(k_sub,:)     = random(pd,ntr,1);
    sub_o(k_sub,:)  = offers(kO(k_sub,:));                                   % resample Offer in pdf (="soft-max")
    
    % compute frequencies
    for k = 0:10
        sub_o_freq(k_sub,k+1)  = sum(sub_o(k_sub,:)==k);
    end
    sub_o_freq(k_sub,:) = sub_o_freq(k_sub,:)./(ntr);
    
    % fit the prior
    n_rep           = 10;
    parameters_rep  = NaN(n_rep,3);
    ll_rep          = NaN(n_rep,1);
    
    
    lb = [-20,0,0]';
    ub = [0,10,10]';
    for k_rep = 1:n_rep
        x0 = [-10*rand(),10*rand(),10*rand()]';
        [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) estimate_priors_cont(x,sub_o(k_sub,:)),x0,[],[],[],[],lb,ub,[],options);
    end
    [~,pos]                  = min(ll_rep);
    parameters(k_sub,:)      = parameters_rep(pos(1),:);
    ll(k_sub)                = ll_rep(pos(1),:);
    
end


figure;
set(gcf,'Color',[1,1,1])
subplot(1,3,1)
hold on
mtp = mean(sub_o_freq);
stp = std(sub_o_freq)./sqrt(nsub);
bar(offers,mtp,'FaceColor',.7.*[1,1,1])
errorbar(offers,mtp,stp,'k','LineStyle','none')
plot(offers,mean(pc),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[-.5 10.5],...
    'XTick',0:10)
xlabel('offers')
ylabel('frequency (%)')
legend('Emprirical','Model')
%
subplot(1,3,2)
plot(offers,mean(EV),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated Expected value')

subplot(1,3,3)
plot(offers,mean(PA),'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')


figure;
set(gcf,'Color',[1,1,1])
for k = 1:3
    subplot(1,3,k)
    plot(PARAMS(:,k),parameters(:,k),'o',...
        'MarkerFaceColor',[1,1,1],...
        'MarkerEdgeColor',[0,0,0])
    
    [R,P] = corr(PARAMS(:,k),parameters(:,k));
end


