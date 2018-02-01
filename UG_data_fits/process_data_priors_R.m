% to make all the matlab formatted data into R format, must be called from
% the directory that the script is saved in
clear
close all
clc

data_path = ['..' filesep 'data' filesep 'processed' ];

priors_soc_dir    = dir([data_path filesep '*_SocPriors*']);
priors_nonsoc_dir = dir([data_path filesep '*NonSocPriors*']);
% % UG_dir            = dir([data_path filesep '*UG*']);

Dict_dir            = dir([data_path filesep '*Dict*']);

cur_st   = 1; 
cur_end  = 96;                % each Dictator has 96 rows (or is supposed to)
all_data = cell(100, 5);      % starting with fewer rows than will be required because I don't know how long it will be

for ss = 1:length(Dict_dir)
    sub_name = Dict_dir(ss).name(1:12);
    flnmSoc = fullfile(data_path,strcat(sub_name,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_path,strcat(sub_name,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
% %         sub_age = load(flnmSoc, 'sub_age');
% %         if isnan(sub_age.sub_age)
% %             sub_age = load(flnmNonSoc, 'sub_age');
% %         end
% %         sub_age = sub_age.sub_age;
        
        % for the priors datafile
        % column 1: offer
        % column 2: option 1
        % column 3: option 2
        % column 4: payment
        % column 5: subject name
        % column 6: trial 
        myfriend = load(flnmSoc, 'subdata');
        
        simple_count = 0;
        trial_count  = 0;
        for hmm = cur_st:cur_st+length(myfriend.subdata)-1
            simple_count     = simple_count + 1;
            trial_count      = trial_count + 1;
            all_data{hmm, 1} = myfriend.subdata(simple_count, 3); % column 3 of Dictator matrix is the offer selected
            all_data{hmm, 2} = myfriend.subdata(simple_count, 1); % option 1
            all_data{hmm, 3} = myfriend.subdata(simple_count, 2); % option 2 
            all_data{hmm, 4} = sub_name;                          % sub_name
            all_data{hmm, 5} = trial_count;                       % trial
            
        end
        
        cur_st = cur_st+length(myfriend.subdata);

    end
    
end

all_header = {'offer', 'option1', 'option2', 'sub_name', 'trial'};
fid = fopen(['..' filesep 'data' filesep 'processed_R' filesep 'Social_priors_all.txt'], 'w' );

for ii = 1:length(all_data)+1 % +1 for header
    if ii == 1
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', all_header{1,:});
    else
        fprintf(fid, '%d\t%d\t%d\t%s\t%d\n', all_data{ii-1, :});
    end
end





% nonsocial priors

cur_st   = 1; 
cur_end  = 96;                % each Dictator has 96 rows (or is supposed to)
all_data = cell(100, 5);      % starting with fewer rows than will be required because I don't know how long it will be

for ss = 1:length(Dict_dir)
    sub_name = Dict_dir(ss).name(1:12);
    flnmSoc = fullfile(data_path,strcat(sub_name,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_path,strcat(sub_name,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
% %         sub_age = load(flnmSoc, 'sub_age');
% %         if isnan(sub_age.sub_age)
% %             sub_age = load(flnmNonSoc, 'sub_age');
% %         end
% %         sub_age = sub_age.sub_age;
        
        % for the priors datafile
        % column 1: offer
        % column 2: option 1
        % column 3: option 2
        % column 4: payment
        % column 5: subject name
        % column 6: trial 

        myfriend = load(flnmNonSoc, 'subdata');
        
        simple_count = 0;
        trial_count  = 0;
        for hmm = cur_st:cur_st+length(myfriend.subdata)-1
            simple_count     = simple_count + 1;
            trial_count      = trial_count + 1;
            all_data{hmm, 1} = myfriend.subdata(simple_count, 3); % column 3 of Dictator matrix is the offer selected
            all_data{hmm, 2} = myfriend.subdata(simple_count, 1); % option 1
            all_data{hmm, 3} = myfriend.subdata(simple_count, 2); % option 2 
            all_data{hmm, 4} = sub_name;                          % sub_name
            all_data{hmm, 5} = trial_count;                       % trial
            
        end
        
        cur_st = cur_st+length(myfriend.subdata);

    end
    
end

all_header = {'offer', 'option1', 'option2', 'sub_name', 'trial'};
fid = fopen(['..' filesep 'data' filesep 'processed_R' filesep 'NonSocial_priors_all.txt'], 'w' );

for ii = 1:length(all_data)+1 % +1 for header
    if ii == 1
        fprintf(fid, '%s\t%s\t%s\t%s\t%s\n', all_header{1,:});
    else
        fprintf(fid, '%d\t%d\t%d\t%s\t%d\n', all_data{ii-1, :});
    end
end


