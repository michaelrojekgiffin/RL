% this script processes the raw otree datafiles and saves each subjects
% data in a readable excel 
clear all
close all
clc

cd ~/Dropbox/RL/data/pilot/3cond_blocked/

raw_dir     = [pwd filesep 'raw'];


% social priors
psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_social*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(psoc_files)
    txtnm = fullfile([raw_dir,filesep,psoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    if sum(strcmp(txt(1, :), 'group.payout_chosen_trial')) > 0
        op1_col       = 24;
        op2_col       = 25;
        offer_col     = 26;
        acc_col       = 31;
        pay_col       = 29;
    else
        op1_col       = 23;
        op2_col       = 24;
        offer_col     = 25;
        acc_col       = 30;
        pay_col       = 28;
    end
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        sub_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), sub_id{sub_counter}));
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
        all_sub_data(sub_counter).subdata  = subdata;
        all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};

        
    end
    
end

% Go through all the data to make sure there's no duplicate subjects in the
% the multiple files that we include, and if there are take the most
% complete version for that subject 
uniquesubcount = 0;
uniqueidx      = unique(sub_id);
for ii = 1:length(uniqueidx)
    repsub         = find(strcmp(sub_id, uniqueidx{ii}));
    
    for kk = 1:length(repsub)
        tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
        all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
    end
    
end

for ii = 1:length(all_subdata_unique)
    sub_id          = uniqueidx{ii};
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('SocPriors_sub',sub_name);
    flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end






%% Non-social priors
raw_dir     = [pwd filesep 'raw'];


% social priors
psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_social*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(psoc_files)
    txtnm = fullfile([raw_dir,filesep,psoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    if sum(strcmp(txt(1, :), 'group.payout_chosen_trial')) > 0
        op1_col       = 24;
        op2_col       = 25;
        offer_col     = 26;
        acc_col       = 31;
        pay_col       = 29;
    else
        op1_col       = 23;
        op2_col       = 24;
        offer_col     = 25;
        acc_col       = 30;
        pay_col       = 28;
    end
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        sub_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), sub_id{sub_counter}));
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
        all_sub_data(sub_counter).subdata  = subdata;
        all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};

        
    end
    
end

% Go through all the data to make sure there's no duplicate subjects in the
% the multiple files that we include, and if there are take the most
% complete version for that subject 
uniquesubcount = 0;
uniqueidx      = unique(sub_id);
for ii = 1:length(uniqueidx)
    repsub         = find(strcmp(sub_id, uniqueidx{ii}));
    
    for kk = 1:length(repsub)
        tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
        all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
    end
    
end

for ii = 1:length(all_subdata_unique)
    sub_id          = uniqueidx{ii};
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('SocPriors_sub',sub_name);
    flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end











%%
% social priors
psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_social*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(psoc_files)
    txtnm = fullfile([raw_dir,filesep,psoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    if sum(strcmp(txt(1, :), 'group.payout_chosen_trial')) > 0
        op1_col       = 24;
        op2_col       = 25;
        offer_col     = 26;
        acc_col       = 31;
        pay_col       = 29;
    else
        op1_col       = 23;
        op2_col       = 24;
        offer_col     = 25;
        acc_col       = 30;
        pay_col       = 28;
    end
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        sub_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), sub_id{sub_counter}));
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
%         flnm = strcat('SocPriors_sub',sub_id{sub_counter});
%         flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
%         
%         save(flnm_to_save,'subdata','sub_id','sub_name')
        
    end
    
end







%% non social priors
pnonsoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_nonsocial*xlsx*']);

sub_name_col  = 2;
for f = 1:length(pnonsoc_files)
    txtnm = fullfile([raw_dir,filesep,pnonsoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    if sum(strcmp(txt(1, :), 'group.payout_chosen_trial')) > 0
        op1_col       = 24;
        op2_col       = 25;
        offer_col     = 26;
        acc_col       = 31;
        pay_col       = 29;
    else
        op1_col       = 23;
        op2_col       = 24;
        offer_col     = 25;
        acc_col       = 30;
        pay_col       = 28;
    end
    
    for subsub = 1:length(all_sub_name)
        sub_id    = all_sub_name{subsub};
        sub_name  = ['00', num2str(subsub)];
        
        subidx    = find(strcmp(txt(:, sub_name_col), sub_id));
        
        subdata   = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
        flnm = strcat('NonSocPriors_sub',sub_name);
        flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
        
        save(flnm_to_save,'subdata','sub_id','sub_name')
        
    end
    
end


%%









cur_dir = pwd;
txtnm = fullfile(cur_dir,'ital.xlsx');

data_dir = fullfile(cur_dir);


[num,txt,raw] = xlsread(txtnm);

subjects = unique(num(:,1));


for k_sub = 1:length(subjects)
    
    i_sub = subjects(k_sub);
    
    sub_lines        = find(num(:,1)==i_sub);
    data             = num(sub_lines,:);
    role             = strcat(txt{sub_lines(1)+1,3});
    
    for k = 1:size(txt,2)
        column_ID{k}        = txt{1,k};
    end
    
    flnm = strcat('DATA_sub',num2str(i_sub));
    flnm_to_save = fullfile(data_dir,flnm);
    
    save(flnm_to_save,'data','role','column_ID')
    
end
