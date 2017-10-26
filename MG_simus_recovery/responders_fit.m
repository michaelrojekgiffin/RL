% this script gets the responses we collected from the responders and
% extracts their parameters, i.e. the slope and intercepts of their logit
% response functions

clear all
close all
clc

offers = 0:20;

cur_dir = pwd;
txtnm = fullfile(cur_dir,'ultimatum_responders.xlsx');
[num,txt,raw] = xlsread(txtnm);

%%
close all
endow_freq = NaN(5, 20);
endow_count = 0;
figure
for ee = 1:5 % because we have 5 different conditions going up to 20 in steps of 5
    endow_lines = find(num(:,6)==endow_count); % get rows for the conditions (specified by starting endowment)
    for oo = offers % get frequency of each offer being accepted
        endow_freq(ee, oo+1)  = sum(num(endow_lines,4) == 1 & num(endow_lines,7) == oo) / sum(num(endow_lines,7) == oo);
    end
    endow_count     =   endow_count+5;
    subplot(5, 1, ee)
    plot(offers, endow_freq(ee, :), '-r', 'Linewidth', 1.5);
    hold on
    bar(endow_freq(ee, :), 'b');
    xlim([offers(1), offers(end)]);
    hold off
end

%% get params
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));


% Compute the p(offers)
fuckPA      = logitp([beta0,beta1],offers);        % compute proba of accepting the offers given current model

n_rep            = 10;
parameters_rep   = NaN(n_rep,3);
ll_rep           = NaN(n_rep,1);

for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:2),ll_rep(k_rep,1)]=fmincon(@(x) logit(x,offers),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0],[Inf Inf],[],options);
end
[~,pos]             = min(ll_rep);
opponent_parameters = parameters_rep(pos(1),:);
ll                  = ll_rep(pos(1),:);


%%
% Get slope and intercept of distribution against which the subject is
% playind (i.e. estimate their true "acceptance function")
opponent_o_freq = zeros(1, 11);
succes_probs = zeros(length(opponent_o_freq), 1); % probability of success
for k = 0:10
    opponent_o_freq(k+1)       = (sum(opponent_o(:)==k)) / length(opponent_o);
    switch predprey
        case 'prey'
            succes_probs(k+1) = sum(opponent_o_freq(1:k+1));
            opponent_role = 'predator';
        case 'predator'                                         % predator only wins if prey invests less
            succes_probs(k+1) = sum(opponent_o_freq(1:k));     % sum all probabilities below chosen offer
            opponent_role = 'prey';
    end
end
empirical_EV = zeros(length(offers), 1);
switch predprey
    case 'prey'
        empirical_EV = (endow - offers).* succes_probs';
    case 'predator'
        empirical_EV = (endow - offers) + ((endow - offers).* succes_probs');
end

% here I'm fitting to the opponent data in order to get the logit
% choice function that subjects were playing against in order to
% use this in the simulation below.
n_rep            = 10;
parameters_rep   = NaN(n_rep,3);
ll_rep           = NaN(n_rep,1);
opponent_o_array = reshape(opponent_o, [], 1);

for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) laplace_priors_priors(x,opponent_o_array),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
end
[~,pos]             = min(ll_rep);
opponent_parameters = parameters_rep(pos(1),:);
ll                  = ll_rep(pos(1),:);





data_dir = fullfile(cur_dir,'data_matlab');
subjects = unique(num(:,4));


for k_sub = 1:length(subjects)
    
    i_sub = subjects(k_sub);
    
    sub_lines        = find(num(:,4)==i_sub);
    data             = num(sub_lines,:);
    role             = strcat(txt{sub_lines(1)+1,2});
    gender           = strcat(txt{sub_lines(1)+1,12});
    sex_night_before = strcat(txt{sub_lines(1)+1,26});
    birth_control    = strcat(txt{sub_lines(1)+1,27});
    food             = strcat(txt{sub_lines(1)+1,28});
    
    for k = 1:size(txt,2)
        column_ID{k}        = txt{1,k};
    end
    
    flnm = strcat('DATA_sub',num2str(i_sub));
    flnm_to_save = fullfile(data_dir,flnm);
    
    save(flnm_to_save,'data','role','gender','sex_night_before','birth_control','food','column_ID')
end



% Get slope and intercept of distribution against which the subject is
% playind (i.e. estimate their true "acceptance function")
opponent_o_freq = zeros(1, 11);
succes_probs = zeros(length(opponent_o_freq), 1); % probability of success
for k = 0:10
    opponent_o_freq(k+1)       = (sum(opponent_o(:)==k)) / length(opponent_o);
    switch predprey
        case 'prey'
            succes_probs(k+1) = sum(opponent_o_freq(1:k+1));
            opponent_role = 'predator';
        case 'predator'                                         % predator only wins if prey invests less
            succes_probs(k+1) = sum(opponent_o_freq(1:k));     % sum all probabilities below chosen offer
            opponent_role = 'prey';
    end
end
empirical_EV = zeros(length(offers), 1);
switch predprey
    case 'prey'
        empirical_EV = (endow - offers).* succes_probs';
    case 'predator'
        empirical_EV = (endow - offers) + ((endow - offers).* succes_probs');
end

% here I'm fitting to the opponent data in order to get the logit
% choice function that subjects were playing against in order to
% use this in the simulation below.
n_rep            = 10;
parameters_rep   = NaN(n_rep,3);
ll_rep           = NaN(n_rep,1);
opponent_o_array = reshape(opponent_o, [], 1);

for k_rep = 1:n_rep
    [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) laplace_priors_priors_MG_2017_10_03(x,opponent_o_array, opponent_role),[10*randn() 10*rand()  10*rand()],[],[],[],[],[-Inf 0 0],[Inf Inf Inf],[],options);
end
[~,pos]             = min(ll_rep);
opponent_parameters = parameters_rep(pos(1),:);
ll                  = ll_rep(pos(1),:);