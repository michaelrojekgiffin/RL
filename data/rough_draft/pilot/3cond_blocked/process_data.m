% this script processes the raw otree datafiles and saves each subjects
% data in a readable excel 
clear all
close all
clc

% cd ~/Dropbox/RL/data/pilot/3cond_interleafed/
cd ~/Dropbox/RL/data/pilot/3cond_blocked/

raw_dir     = [pwd filesep 'raw'];


%% social priors
psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_social*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(psoc_files)
    txtnm = fullfile([raw_dir,filesep,psoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    op1_col       = find(strcmp(txt(1, :), 'group.option1'));
    op2_col       = find(strcmp(txt(1, :), 'group.option2'));
    offer_col     = find(strcmp(txt(1, :), 'group.amount_offered'));
    acc_col       = find(strcmp(txt(1, :), 'group.offer_accepted'));
    pay_col       = find(strcmp(txt(1, :), 'group.payout'));
    
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
        proper_id{uniquesubcount} = uniqueidx{ii};
    end
    
end

for ii = 1:length(all_subdata_unique)
    sub_id          = proper_id{ii};
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('Sub_', sub_id, '_SocPriors');
    flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end


%% Non-social priors
clear all
raw_dir     = [pwd filesep 'raw'];

% social priors
pnonsoc_files    = dir([raw_dir,filesep, '*my_ultimatum_priors_nonsocial*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(pnonsoc_files)
    txtnm = fullfile([raw_dir,filesep,pnonsoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    op1_col       = find(strcmp(txt(1, :), 'group.option1'));
    op2_col       = find(strcmp(txt(1, :), 'group.option2'));
    offer_col     = find(strcmp(txt(1, :), 'group.amount_offered'));
    acc_col       = find(strcmp(txt(1, :), 'group.offer_accepted'));
    pay_col       = find(strcmp(txt(1, :), 'group.payout'));
    
    
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
        proper_id{uniquesubcount} = uniqueidx{ii};
    else
        fprintf('%s is an empty matrix, row %d of all_sub_data.\n', all_sub_data(repsub(kk)).sub_id, repsub(kk));
    end
    
end


for ii = 1:length(all_subdata_unique)
    sub_id          = proper_id{ii};
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('Sub_', sub_id, '_NonSocPriors');
    flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end




%% UG
clear all
raw_dir     = [pwd filesep 'raw'];

% psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_3cond_interleafed*xlsx*']);
psoc_files    = dir([raw_dir,filesep, '*my_ultimatum_3cond_block*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(psoc_files)
    txtnm         = fullfile([raw_dir,filesep,psoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name  = unique(txt(2:end, sub_name_col)); 
    
    offer_col     = find(strcmp(txt(1, :), 'group.amount_offered'));
    acc_col       = find(strcmp(txt(1, :), 'group.offer_accepted'));
    pay_col       = find(strcmp(txt(1, :), 'group.payout'));
    
    social_col    = find(strcmp(txt(1, :), 'group.soc_or_no'));        % string
    opponent_col  = find(strcmp(txt(1, :), 'group.current_opponent')); % string
    
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        sub_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), sub_id{sub_counter}));
        
        soc_array = NaN(length(subidx), 1);
        opp_array = NaN(length(subidx), 1);
        for cc = 1:length(txt(subidx))
            if strcmp(txt(subidx(cc), social_col), 'Social')
                soc_array(cc) = 1;
            else
                soc_array(cc) = 0;
            end
            if strcmp(txt(subidx(cc), opponent_col), 'triangle')
                opp_array(cc) = 1;
            elseif strcmp(txt(subidx(cc), opponent_col), 'circle') 
                opp_array(cc) = 2;
            elseif strcmp(txt(subidx(cc), opponent_col), 'square')
                opp_array(cc) = 3;
            end
            
            
        end
        
        subdata                  = [num(subidx-1, [offer_col, acc_col, pay_col]), soc_array, opp_array]; % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
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
        proper_id{uniquesubcount} = uniqueidx{ii};
    end
    
end


for ii = 1:length(all_subdata_unique)
    sub_id          = proper_id{ii};
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('Sub_', sub_id, '_UG');
    flnm_to_save = fullfile([pwd, filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end


