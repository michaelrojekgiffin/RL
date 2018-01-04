% to make all the matlab formatted data into R format
clear
close all
clc

data_path = ['..' filesep 'data' filesep 'processed' ];

priors_soc_dir    = dir([data_path filesep '*_SocPriors*']);
priors_nonsoc_dir = dir([data_path filesep '*NonSocPriors*']);
UG_dir            = dir([data_path filesep '*UG*']);

cur_st   = 1; 
cur_end  = 240; % each UG has 240 rows
all_data = cell(24000, 7);

for ss = 1:length(UG_dir)
    sub_name = UG_dir(ss).name(1:12);
    flnmSoc = fullfile(data_path,strcat(sub_name,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_path,strcat(sub_name,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        sub_age = load(flnmSoc, 'sub_age');
        if isnan(sub_age.sub_age)
            sub_age = load(flnmNonSoc, 'sub_age');
        end
        sub_age = sub_age.sub_age;
        
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
        myfriend = load(fullfile(data_path, UG_dir(ss).name), 'subdata');
        
        simple_count = 0;
        trial_count  = 0;
        for hmm = cur_st:cur_end
            simple_count     = simple_count + 1;
            trial_count      = trial_count + 1;
            all_data{hmm, 1} = myfriend.subdata(simple_count, 1);
            all_data{hmm, 2} = myfriend.subdata(simple_count, 2);
            all_data{hmm, 3} = myfriend.subdata(simple_count, 3);
            all_data{hmm, 4} = myfriend.subdata(simple_count, 4);
            all_data{hmm, 5} = myfriend.subdata(simple_count, 5);
            all_data{hmm, 6} = myfriend.subdata(simple_count, 6);
            all_data{hmm, 7} = sub_name;
            all_data{hmm, 8} = trial_count;
            all_data{hmm, 9} = sub_age;
            
            if trial_count == 120
                trial_count = 0;
            end
        end
        
        cur_st                        = cur_st+240;
        cur_end                       = cur_end+240;
    end
    
end

all_header = {'offer', 'accepted', 'payoff', 'social', 'opponent', 'explicit', 'sub_name', 'trial', 'age'};
fid = fopen(['..' filesep 'data' filesep 'processed_R' filesep 'UG_all.txt'], 'w' );

for ii = 1:length(all_data)+1 % +1 for header
    if ii == 1
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', all_header{1,:});
    else
        fprintf(fid, '%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n', all_data{ii-1, :});
    end
end
