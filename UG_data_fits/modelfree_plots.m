% to make some model free plots corroborating the statistics I've found in
% R. The first part of this script is actually reshaping the data, and it
% would be a good idea to add this to the process_data_R script to just
% save it in a better format

clear
close all force
clc


% standard error function
ste = @(x) (std(x))/(sqrt(length(x)));
% ste = @(x) (std(x))/(sqrt(length(x)));

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% find paths and import data
% this gets all subjects into one big matrix, not super elegant but closer
% to what I'm used to in R so should help speed things along
cur_dir = pwd;
%  project_name = 'RL_PreyPredator';
project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);

data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed_R', 'UG.mat');

load(data_dir);
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% for the UG datafile
% column 1: offer
% column 2: accepted (1 for yes, 0 for no)
% column 3: payment
% column 4: social (1 for social, 0 for non-social)
% column 5: opponent, 1 for starting endowment of 0, 2 for starting endowment of 10, and 3 for starting endowment of 20
% column 6: explicit (1 for explicit, 0 for implicit)
% column 7: subject name
% column 8: trial
% column 9: subject age
all_header = {'offer', 'accepted', 'payoff', 'social', 'opponent', 'explicit', 'sub_name', 'trial', 'age'};
%% format data in nice way
ntr = 120; % 25 is the number of trials for which all subjects are exposed to each opponent

soc     = UG(UG(:,4) == 1, :);

soc_imp = soc(soc(:,6) == 0, :);
soc_exp = soc(soc(:,6) == 1, :);

nonsoc  = UG(UG(:,4) == 0, :);

nonsoc_exp = nonsoc(nonsoc(:,6) == 1, :);
nonsoc_imp = nonsoc(nonsoc(:,6) == 0,  :);

explicit = UG(UG(:,6) == 1, :);
implicit = UG(UG(:,6) == 0, :);

% for ugm I want a nXmXp matrix where n is condtion (of which we have 4),
% m is opponent (of which we have 3), and p is trial (of which we have a
% variable number, but I will cap it at 40 because that should be the
% average sincel that's the number of total trails divided by the numner of
% opponents)
all_subs_i = unique(implicit(:, 7));
all_subs_e = unique(explicit(:, 7));
nsub = length(all_subs_i)+length(all_subs_e);

% order of dimensions are:
% 1: subject
% 2: social nonsocial
% 3: opponent
% 4: trial
ugm = NaN(length(all_subs_i), 4, 3, ntr);
ugs = NaN(length(all_subs_i), 4, 3, ntr);

% ugs = NaN(4, 3, ntr, length(unique(nonsoc(:, 7))));

for ss = 1:nsub/2 % nsub/2 because half subjects are in explicit and half are implicit
    for oo = 1:size(ugm, 3)
        
        tmp1 = (soc_imp(soc_imp(:, 5) == oo & soc_imp(:, 7) == all_subs_i(ss), 1));
        tmp2 = (nonsoc_imp(nonsoc_imp(:, 5) == oo & nonsoc_imp(:, 7) == all_subs_i(ss), 1));
        
        tmp3 = (soc_exp(soc_exp(:, 5) == oo & soc_exp(:, 7) == all_subs_e(ss), 1));
        tmp4 = (nonsoc_exp(nonsoc_exp(:, 5) == oo & nonsoc_exp(:, 7) == all_subs_e(ss), 1));
       
        ugm(ss, 1, oo, 1:length(tmp1)) = tmp1;
        ugm(ss, 2, oo, 1:length(tmp2)) = tmp2;
        ugm(ss, 3, oo, 1:length(tmp3)) = tmp3;
        ugm(ss, 4, oo, 1:length(tmp4)) = tmp4;
        
    end
end

%% plots that work!
close all

figure
plot_titles = {'SocImp', 'NonSocImp', 'SocExp', 'NonSocExp'};
plot_leg = {'0 endow', '10 endow', '20 endow'};
all_leg = {'r-^', 'b-.o', 'k:s'};


for pp = 1:4 
   
    subplot(2, 2, pp)
    
    hold on
    for ii = 1:3
        
        %     plot(1:ntr, mean(squeeze((ugm(:, 2, ii, :)))), all_leg{ii}, 'markers',5, 'linewidth', 1)
        
%         errorbar(mean(squeeze((ugm(:, pp, ii, :)))), ste(squeeze((ugm(:, pp, ii, :)))), all_leg{ii}, 'markers',8, 'linewidth', 1)
        
        ugm_tmp = squeeze((ugm(:, pp, ii, :))); % subject X trial matrix
        
        umeans = [];
        ustes  = [];
        for nn = 1:size(ugm_tmp, 2)
            ugm_clean = reshape(ugm_tmp(:,nn), [], 1);
            
            ugm_clean(any(isnan(ugm_clean), 2), :) = [];
            
            if ~isempty(ugm_clean)
                umeans(nn) = mean(ugm_clean);
                ustes(nn)  = ste(ugm_clean);
            end
        end
       
        errorbar(umeans, ustes, all_leg{ii}, 'markers',8, 'linewidth', 1)
        
    end
    ylim([0 12])
    xlim([0 50])
    legend(plot_leg, 'location', 'northwest')
    title(plot_titles{pp})
    hold off
end
   