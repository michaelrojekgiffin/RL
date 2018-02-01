% to make all the matlab formatted data into R format, must be called from
% the directory that the script is saved in
clear
close all
clc

data_path = ['..' filesep 'data' filesep 'processed' ];

priors_soc_dir    = dir([data_path filesep '*_SocPriors*']);
priors_nonsoc_dir = dir([data_path filesep '*NonSocPriors*']);
% % UG_dir            = dir([data_path filesep '*UG*']);

% Dict_dir            = dir([data_path filesep '*Dict*']);
Risk_dir            = dir([data_path filesep '*Risk*']);

cur_st   = 1; 
cur_end  = 96;                % each Dictator has 96 rows (or is supposed to)
all_data = cell(98, 2);      % starting with fewer rows than will be required because I don't know how long it will be

for ss = 1:length(Risk_dir)
    sub_name = Risk_dir(ss).name(1:12);
    flnmSoc = fullfile(data_path,strcat(sub_name,'_SocPriors.mat'));
    flnmNonSoc = fullfile(data_path,strcat(sub_name,'_NonSocPriors.mat'));
    
    if exist(flnmSoc,'file') && exist(flnmNonSoc,'file')
        
        % for the priors datafile
        % column 1: risk value
        % column 2: subject name

        myfriend = load(fullfile(data_path, Risk_dir(ss).name), 'subdata');
        
        simple_count = 0;
        trial_count  = 0;
        for hmm = cur_st:cur_st+length(myfriend.subdata)-1
            simple_count     = simple_count + 1;
            trial_count      = trial_count + 1;
            all_data{hmm, 1} = myfriend.subdata(simple_count, 1); 
            all_data{hmm, 2} = sub_name;                          % sub_name
            
        end
        
        cur_st = cur_st+length(myfriend.subdata);

    end
    
end

all_header = {'risk', 'sub_name'};
fid = fopen(['..' filesep 'data' filesep 'processed_R' filesep 'risk_all.txt'], 'w' );

for ii = 1:length(all_data)+1 % +1 for header
    if ii == 1
        fprintf(fid, '%s\t%s\n', all_header{1,:});
    else
        fprintf(fid, '%d\t%s\n', all_data{ii-1, :});
    end
end

