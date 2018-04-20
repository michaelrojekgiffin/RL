% performs model comparison on first tier of model fitting (first 4 models)

%% load priors
load('big_hor_subject_fit_workspace.mat')


% importantly, the variable LEARN_LPP has the dimensions
% 1: subject
% 2: model
% 3: role (1 is predator, 2 is prey) - must remove 0 values before
% comparing models
%% find paths
cur_dir = pwd;
% project_name = 'RL_PreyPredator';
project_name = 'RL/PPG'; % for use in michael's dropbox

study_name   = 'matlab';

findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,['data_',study_name]);

fl_list = dir(strcat(data_dir,filesep,'*DATA*.mat*'));
nsub = length(fl_list);


%% some params
nfpm = [2,3,2,3];
Nmodel = 2;
ntr = 60;
ncond = 1;

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:10;
endow  = 10*ones(1,numel(offers));% parameters of the simulation

%% model comparison
options.families = {[1,2], [3,4]} ;

bic = [];
X = pred_LPP;
for k_mod= 1:4
    bic(:,k_mod)=-2*-X(:,k_mod) + nfpm(k_mod)*log(60); % l2 is already positive
end
[postBMC1,outBMC1]=VBA_groupBMC(-bic'./2);

bic = [];
X = prey_LPP;
for k_mod= 1:4
    bic(:,k_mod)=-2*-X(:,k_mod) + nfpm(k_mod)*log(60); % l2 is already positive
end
[postBMC1,outBMC1]=VBA_groupBMC(-bic'./2);




