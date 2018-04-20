% performs model comparison on first tier of model fitting (first 4 models)

%% load priors
load('Learning_modelspace1_2017_09_03.mat')


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
nfpm = [4,5,4,5];
Nmodel = 2;
ntr = 60;
ncond = 1;

%% set params and functions
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
offers = 0:1:10;
endow  = 10*ones(1,numel(offers));% parameters of the simulation

%% model comparison
options.families = {[1,2], [3,4]} ;

pred_LPP = squeeze(LEARN_LPP(:,:,1));
pred_LPP(~any((pred_LPP), 2), :) = [];

[postBMC1,outBMC1]=VBA_groupBMC(-pred_LPP');


prey_LPP = squeeze(LEARN_LPP(:,:,2));
prey_LPP(~any((prey_LPP), 2), :) = [];

[postBMC1,outBMC1]=VBA_groupBMC(-prey_LPP');

% models 2 and 4 are clear winners but indistinguishable from one-another

all_LPP = [pred_LPP; prey_LPP];
[postBMC1,outBMC1]=VBA_groupBMC(-all_LPP');


%% if we punish the extra free parameters again with bic, then 1 and 3 are tied for the win
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



%% Let's just go with model 2 for now...
% I guess we'll need a more sensitive way to see if models 2 and 4 actually
% differ from one-another. First I should run some simulations to see if
% they're identifiable with the PPG
pred_params = squeeze(LEARN_parametersLPP(:,:,3,2));
pred_params(~any((pred_params), 2), :) = [];

histogram(pred_params(:, 2))
% all the learning rates are .5....




