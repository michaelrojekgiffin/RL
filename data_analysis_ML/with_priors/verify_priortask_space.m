clear
close all force
clc

load('Priors_test1_2017_12_01')

%% set dir
cur_dir = pwd;
project_name = 'RL_PreyPredator';
% project_name = 'RL'; % for use in michael's dropbox

%% localize data
findnm = strfind(pwd,project_name);
data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','implicit','processed');

fl_list = dir(strcat(data_dir,filesep,'Sub_*_SocPriors.mat'));

%% load and set params

offers  = 0:1:20;

good_sub = NaN(nsub,1);

%% subject loop
for k_sub = 1:nsub;
    
    sub_nm = fl_list(k_sub).name(5:12);
    flnmSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_dir,strcat('Sub_',sub_nm,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        good_sub(k_sub) = 1;
        
        for kSoc = 1:2
            %% Get UG data
            
            switch kSoc
                case 1
                    load(flnmNonSoc)
                case 2
                    load(flnmSoc)
            end
            
       chO = sortrows(subdata(:,[1,2]));
       
       chO_mat(k_sub,:,:) = chO;
            
        end
    end
    
end