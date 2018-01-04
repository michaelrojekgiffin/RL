% this script counts how many subjects are in the conditno
clear
close all
clc

%------------------------------------------
% change depending on where the subject's data is stored
%------------------------------------------
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '3cond_interleafed' filesep 'processed' ];
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '3cond_blocked' filesep 'processed' ];
data_path = ['..' filesep 'data' filesep 'processed' ];
% data_path = ['..' filesep 'data' filesep 'pilot' filesep '60_trials' filesep '3cond_blocked' filesep 'processed' ];
%------------------------------------------
%------------------------------------------

priors_soc_dir    = dir([data_path filesep '*_SocPriors*']);
priors_nonsoc_dir = dir([data_path filesep '*NonSocPriors*']);
UG_dir            = dir([data_path filesep '*UG*']);

all_there_count = 0;
for ss = 1:length(UG_dir)
    sub_name = UG_dir(ss).name(1:12);
    
    all_there = 0;
    for fsub = 1:length(priors_soc_dir)
        tmpidx = (strcmp(priors_soc_dir(fsub).name(1:12), sub_name));
        if tmpidx == 1
            soc_priors_idx = tmpidx;
            all_there = all_there+1;
        end
    end
    for fsub = 1:length(priors_nonsoc_dir)
        tmpidx = (strcmp(priors_nonsoc_dir(fsub).name(1:12), sub_name));
        if tmpidx == 1
            nonsoc_priors_idx = tmpidx;
            all_there = all_there+1;
        end
    end
    if all_there == 2
        all_there_count = all_there_count + 1;
    end
    % only run the following if the subject has both priors and UG data
end
all_there_count