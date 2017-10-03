% what I need to do here is, first, get the first investments of the entire
% dataset (can do it crudely in R), and then get the priors based on
% everything, save them to matlab files (this script already does it, make
% sure I replace the functions here with the Laplace estimation scripts),
% then I can fit the models to my actual data!
% clc
clear
close all

cur_dir     = pwd;
data_dir    = fullfile(cur_dir,'data_matlab');
fl_dir      = dir(strcat(data_dir,filesep,'DATA_sub*'));
nsub        = length(fl_dir);
% sub_o       = NaN(nsub, 1)./2; % preallocation needs to account for it
% being just prey

ntr = 1; % number of trials used per subject(1:first, 2:first two, 3:first three, etc...)

%% for prey
% I got the below distribution from the first investments made by subjects
% playing as prey, in priors_dist.R in 
% /Carsten PhD/hormones/scripts/modeling
prey_prior_dist = [9  6  5  4  8  8 10  5  7  6  5  8  6  5  7  5  5  6  7  8  4  6  9  7  5  5  5  9  2  6  7  6  6  5  8  4  6  7  8  9  9  7  9  8  9  5  5  8  4  5  6  3  6  6  4  8 2  7  9  7  8 10  8  4  8 10  7  6  7  5  8  8  5  5  8  6  6  6  8  7  7  8  7  5  6  4  6 10  3  0  1  4  6  7  2  7  6  5  5  4  7  8  6  8  8  8  6  5  0  4  8  6 8  6  7  5  5  6  5  7  0  6  6  8  5  6  5  4  5  4];

k_prey = 0;
sub_o = [];
for k_sub   = 1:nsub
    
    flnm    = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
%             k_prey                  = k_prey+1;
%             sub_o(k_prey, 1:ntr)    = data(1:ntr,5);
        case 'prey'       
            k_prey                  = k_prey+1;
            sub_o(k_prey, 1:ntr)    = data(1:ntr,5);    
    end
end

predprey = 'prey';

% compute frequences
sub_o = sub_o(:);
sub_o_freq = zeros(1, 11);
for k = 0:10
    sub_o_freq(k+1)  = sum(sub_o==k);
end
sub_o_freq = sub_o_freq./(k_prey.*ntr);

% fit the prior
options         = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);
n_rep           = 10;
parameters_rep  = NaN(n_rep,3);
ll_rep          = NaN(n_rep,1);
for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) PriorEstimation_mg_2017_09_04(x,sub_o, predprey),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
end
[~,pos]         = min(ll_rep);
parameters      = parameters_rep(pos(1),:);
ll              = ll_rep(pos(1),:);

% Save parameters .mat file
flnm = 'prey_priors';
flnm_to_save = fullfile(data_dir,flnm);
save(flnm_to_save, 'parameters')

% Plot Everything
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

PA      = logitp([parameters(1),parameters(2)],offers);            % compute proba of accepting the offers given current model
EV      = (endow - offers).* PA;                                   % compute EV of the offers given current model
pc = exp(parameters(3).*EV) ./ sum(exp(parameters(3).*EV));


figure;
set(gcf,'Color',[1,1,1])
subplot(1,3,1)
title('Prey')
hold on
bar(offers,sub_o_freq,'FaceColor',.7.*[1,1,1])
plot(offers,pc,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[-.5 10.5],...
    'XTick',0:10)
xlabel('offers')
ylabel('frequency (%)')
legend('Emprirical','Model')

subplot(1,3,2)
plot(offers,EV,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated Expected value')


subplot(1,3,3)
plot(offers,PA,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')

%% for predator
k_pred = 0;
sub_o = [];
for k_sub   = 1:nsub
    
    flnm    = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
            k_pred                  = k_pred+1;
            sub_o(k_pred, 1:ntr)    = data(1:ntr,5);  
        case 'prey'       
%             k_pred                  = k_pred+1;
%             sub_o(k_pred, 1:ntr)    = data(1:ntr,5);    
    end
end
predprey = 'predator';
% compute frequences
sub_o = sub_o(:);
sub_o_freq = zeros(1, 11);
for k = 0:10
    sub_o_freq(k+1)  = sum(sub_o==k);
end
sub_o_freq = sub_o_freq./(k_pred.*ntr);

% fit the prior
options         = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);
n_rep           = 10;
parameters_rep  = NaN(n_rep,3);
ll_rep          = NaN(n_rep,1);
for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) PriorEstimation_mg_2017_09_04(x,sub_o, predprey),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
end
[~,pos]         = min(ll_rep);
parameters      = parameters_rep(pos(1),:);
ll              = ll_rep(pos(1),:);

% Save parameters .mat file
flnm = 'predator_priors';
flnm_to_save = fullfile(data_dir,flnm);
save(flnm_to_save, 'parameters')


% Plot Everything
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

PA      = logitp([parameters(1),parameters(2)],offers);            % compute proba of accepting the offers given current model
EV      = (endow - offers) + ((endow - offers).* PA);              % compute EV of the offers given current model
pc = exp(parameters(3).*EV) ./ sum(exp(parameters(3).*EV));


figure;
set(gcf,'Color',[1,1,1])
subplot(1,3,1)
title('Predator')
hold on
bar(offers,sub_o_freq,'FaceColor',.7.*[1,1,1])
plot(offers,pc,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
set(gca,'XLim',[-.5 10.5],...
    'XTick',0:10)
xlabel('offers')
ylabel('frequency (%)')
legend('Emprirical','Model')

subplot(1,3,2)
plot(offers,EV,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated Expected value')


subplot(1,3,3)
plot(offers,PA,'-o',...
    'Color',0.*[1,1,1],...
    'MarkerFaceColor',1.*[1,1,1],...
    'MarkerEdgeColor',0.*[1,1,1])
xlabel('offers')
ylabel('Estimated probability of Acceptance')