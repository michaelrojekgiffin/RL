clc
clear
close all

cur_dir     = pwd;
data_dir    = fullfile(cur_dir,'data_matlab');
fl_dir      = dir(strcat(data_dir,filesep,'DATA_sub*'));
nsub        = length(fl_dir);
% sub_o       = NaN(nsub, 1)./2; % preallocation needs to account for it
% being just prey

ntr = [1 30 60]; % nb of trials used per subject(1:first, 2:first two, 3:first three, etc...)
ntest = length(ntr);

k_prey = 0;
for k_sub   = 1:nsub
    
    flnm    = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
        case 'prey'
            k_prey                  = k_prey+1;
            sub_o(k_prey, 1:ntest)    = data(ntr,5);
    end
end


% fit the prior
options         = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);
n_rep           = 10;
parameters_rep  = NaN(n_rep,7);
ll_rep          = NaN(n_rep,1);
for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:7),ll_rep(k_rep,1)]=fmincon(@(x) LearnEstimation(x,sub_o),[5*rand(1,7)],[],[],[],[],[repmat([-Inf 0],1,3) 0],[repmat(Inf,1,7)],[],options);
end
[~,pos]       = min(ll_rep);
parameters    =   parameters_rep(pos(1),:);
ll            =   ll_rep(pos(1),:);



figure(1);
set(gcf,'Color',[1,1,1])

figure(2);
set(gcf,'Color',[1,1,1])

sub_o_freq = zeros(ntest, 11);
for k_test = 1:ntest
    
    a = k_test./ntest;
    b = (ntest - k_test)./ntest;
    
    for k = 0:10
        sub_o_freq(k_test,k+1)  = sum(sub_o(:,k_test)==k)./k_prey;
    end
    
    % Plot Everything
    offers  = 0:1:10;
    endow   = 10*ones(1,numel(offers));
    logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
    
    PA      = logitp([parameters((k_test-1)*2+1),parameters(k_test*2)],offers); % compute proba of accepting the offers given current model
    EV      = (endow - offers).* PA;                                            % compute EV of the offers given current model
    pc = exp(parameters(7).*EV) ./ sum(exp(parameters(7).*EV));
    
    
    
    figure(1)
    subplot(1,3,k_test)
    hold on
    bar(offers,sub_o_freq(k_test,:),'FaceColor',.7.*[1,1,1])
    plot(offers,pc,'-o',...
        'Color',0.*[1,1,1],...
        'MarkerFaceColor',1.*[1,1,1],...
        'MarkerEdgeColor',0.*[1,1,1])
    set(gca,'XLim',[-.5 10.5],...
        'XTick',0:10)
    xlabel('offers')
    ylabel('frequency (%)')
    legend('Emprirical','Model')
    
    
    figure(2)
    subplot(1,2,1)
    hold on
    plot(offers,EV,'-o',...
        'Color',[a,0,b],...
        'MarkerFaceColor',[a,0,b],...
        'MarkerEdgeColor',[a,0,b])
    
    
    
    subplot(1,2,2)
    hold on
    plot(offers,PA,'-o',...
        'Color',[a,0,b],...
        'MarkerFaceColor',[a,0,b],...
        'MarkerEdgeColor',[a,0,b])
    
    
end


figure(2)
subplot(1,2,1)
xlabel('offers')
ylabel('Estimated Expected value')
legend('Trial 1','Trial 30','Trial 60')

subplot(1,2,2)
xlabel('offers')
ylabel('Estimated probability of Acceptance')
legend('Trial 1','Trial 30','Trial 60')
