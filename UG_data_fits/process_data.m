% this script processes the raw otree datafiles and saves each subjects
% data in a readable excel 
clear all
close all
clc

raw_dir     = ['..', filesep, 'data' filesep,'raw'];


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
    age_col       = find(strcmp(txt(1, :), 'player.age'));
    
    test_ques_imp = find(strcmp(txt(1, :), 'player.test_question_8')); % this question exists for implicit and explicit
    test_ques_exp = find(strcmp(txt(1, :), 'player.test_question_9')); % this questions only exists for explicit
    % if test_ques_imp is not empty and test_ques_exp is empty, it's
    % implicit, if neither are empty it's the explicit condition, and if
    % both are empty it means that non-social priors included the
    % instructions (it was couter-balanced)
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
%         sub_id{sub_counter}      = all_sub_name{subsub};
        tmp_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), tmp_id{sub_counter}));
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        sub_age                  = num(subidx(1)-1, age_col); 
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
        % check for implicit or explicit conditions
        if isempty(txt(subidx(1), test_ques_exp))
            sub_id{sub_counter}      = all_sub_name{subsub};
            all_sub_data(sub_counter).subdata  = subdata;
            all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
            all_sub_data(sub_counter).sub_age  = sub_age;
            all_sub_data(sub_counter).explicit = false;
        else 
            if strcmp(txt(subidx(1), test_ques_imp), 'NON-SOCIAL') && strcmp(txt(subidx(1), test_ques_exp), '')
%                 fprintf('%s non-social filled exp empty \n', all_sub_name{subsub});
                sub_id{sub_counter}      = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = false;
            elseif strcmp(txt(subidx(1), test_ques_imp), '') && strcmp(txt(subidx(1), test_ques_exp), '')
%                 fprintf('%s both empty social priors \n', all_sub_name{subsub});
                sub_id{sub_counter}      = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = false;
            else
                %  fprintf('%s is in explicit condition social priors \n', all_sub_name{subsub});
                sub_id{sub_counter}      = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = true;
            end
        end
    end
    
end

% Go through all the data to make sure there's no duplicate subjects in the
% the multiple files that we include, and if there are take the most
% complete version for that subject 
uniquesubcount = 0;
uniqueidx      = unique(sub_id);
for ii = 1:length(uniqueidx)
    repsub         = find(strcmp(sub_id, uniqueidx{ii}));
    
%     for kk = 1:length(repsub)
%         tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
%         all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
%         tmp(kk).sub_age = all_sub_data(repsub(kk)).sub_age;
%         tmp(kk).explicit = all_sub_data(repsub(kk)).explicit;
%     end
%     
%     % then take the biggest one
%     if ~isempty(all_sub_data(repsub(kk)).subdata)
%         uniquesubcount = uniquesubcount +1;
%         bigidx = find(max(all_lengths));
%         all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
%         all_subdata_unique(uniquesubcount).sub_age = tmp(bigidx).sub_age;
%         all_subdata_unique(uniquesubcount).explicit = tmp(bigidx).explicit;
%     end
    
    
    for kk = 1:length(repsub)
        tmp(kk).data     = all_sub_data(repsub(kk)).subdata;
        all_lengths(kk)  = length(all_sub_data(repsub(kk)).subdata);
        tmp(kk).sub_age  = all_sub_data(repsub(kk)).sub_age;
        tmp(kk).explicit = all_sub_data(repsub(kk)).explicit;
        tmp(kk).sub_id  = uniqueidx{ii};
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata  = tmp(bigidx).data;
        all_subdata_unique(uniquesubcount).sub_age  = tmp(bigidx).sub_age;
        all_subdata_unique(uniquesubcount).explicit = tmp(bigidx).explicit;
        all_subdata_unique(uniquesubcount).sub_id = tmp(bigidx).sub_id;
    else
        fprintf('%s is an empty matrix, row %d of all_sub_data.\n', all_sub_data(repsub(kk)).sub_id, repsub(kk));
    end
    
    
end

for ii = 1:length(all_subdata_unique)
%     sub_id          = uniqueidx{ii};
    sub_id          = all_subdata_unique(ii).sub_id;
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    sub_age         = all_subdata_unique(ii).sub_age;
    explicit        = all_subdata_unique(ii).explicit;
    
    flnm            = strcat('Sub_', sub_id, '_SocPriors');
    flnm_to_save = fullfile(['..', filesep, 'data', filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name', 'sub_age', 'explicit')
end


%% Non-social priors
clear all
raw_dir     = ['..', filesep, 'data' filesep,'raw'];

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
    age_col       = find(strcmp(txt(1, :), 'player.age'));
    test_ques_imp = find(strcmp(txt(1, :), 'player.test_question_8')); % this question exists for implicit and explicit
    test_ques_exp = find(strcmp(txt(1, :), 'player.test_question_9')); % this questions only exists for explicit
    % if test_ques_imp is not empty and test_ques_exp is empty, it's
    % implicit, if neither are empty it's the explicit condition, and if
    % both are empty it means that social priors included the
    % instructions (it was couter-balanced)
    
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        %         sub_id{sub_counter}      = all_sub_name{subsub};
        tmp_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), tmp_id{sub_counter}));
        sub_age                  = num(subidx(1)-1, age_col); 
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, acc_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
         % this makes sure that I'm only adding the implicit condition
         if isempty(txt(subidx(1), test_ques_exp))
             sub_id{sub_counter}      = all_sub_name{subsub};
             all_sub_data(sub_counter).subdata  = subdata;
             all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
             all_sub_data(sub_counter).sub_age  = sub_age;
             all_sub_data(sub_counter).explicit = false;
         else
            if strcmp(txt(subidx(1), test_ques_imp), 'NON-SOCIAL') && strcmp(txt(subidx(1), test_ques_exp), '')
%                 fprintf('%s nonsocial filled exp empty nonsocial priors implicit \n', all_sub_name{subsub});
                sub_id{sub_counter}      = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = false;
            elseif strcmp(txt(subidx(1), test_ques_imp), '') && strcmp(txt(subidx(1), test_ques_exp), '')
%                 fprintf('%s both empty nonsocial priors \n', all_sub_name{subsub});
                sub_id{sub_counter}      = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = false;
            else
                %                 fprintf('%s is in explicit condition nonsocial priors \n', all_sub_name{subsub});
                %                 sub_counter = sub_counter-1;
                sub_id{sub_counter}         = all_sub_name{subsub};
                all_sub_data(sub_counter).subdata  = subdata;
                all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};
                all_sub_data(sub_counter).sub_age  = sub_age;
                all_sub_data(sub_counter).explicit = true;
            end
         end
%         all_sub_data(sub_counter).subdata  = subdata;
%         all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};

        
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
        tmp(kk).data     = all_sub_data(repsub(kk)).subdata;
        all_lengths(kk)  = length(all_sub_data(repsub(kk)).subdata);
        tmp(kk).sub_age  = all_sub_data(repsub(kk)).sub_age;
        tmp(kk).explicit = all_sub_data(repsub(kk)).explicit;
        tmp(kk).sub_id  = uniqueidx{ii};
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata  = tmp(bigidx).data;
        all_subdata_unique(uniquesubcount).sub_age  = tmp(bigidx).sub_age;
        all_subdata_unique(uniquesubcount).explicit = tmp(bigidx).explicit;
        all_subdata_unique(uniquesubcount).sub_id = tmp(bigidx).sub_id;
    else
        fprintf('%s is an empty matrix, row %d of all_sub_data.\n', all_sub_data(repsub(kk)).sub_id, repsub(kk));
    end
    
     
    
end


for ii = 1:length(all_subdata_unique)
%     sub_id          = uniqueidx{ii};
    sub_id          = all_subdata_unique(ii).sub_id;
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    sub_age         = all_subdata_unique(ii).sub_age;
    explicit        = all_subdata_unique(ii).explicit;
    
    flnm            = strcat('Sub_', sub_id, '_NonSocPriors');
    flnm_to_save = fullfile(['..', filesep, 'data', filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name', 'sub_age', 'explicit')
end




%% UG
clear all
raw_dir     = ['..', filesep, 'data' filesep,'raw'];

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
        
        % this is needed to deterine if it's in the explicit or implicit
        % condition. 
        
        
        try
            explitive(1)        = load(fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_NonSocPriors.mat']), 'explicit');
            explitive(2)        = load(fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_SocPriors.mat']), 'explicit');
            explicit            = sum([explitive(1).explicit, explitive(2).explicit]);
            
        catch
            if ~exist(fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_NonSocPriors.mat'])) || ~exist(fullfile('..', filesep, 'data', 'processed', filesep, ['Sub_' sub_id{sub_counter} '_SocPriors.mat']))
                fprintf('%s does not exist\n', fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_NonSocPriors.mat']));
                explicit = -999;
            else
                explitive(1)        = load(fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_NonSocPriors.mat']), 'explicit');
                explitive(2)        = load(fullfile('..', filesep, 'data', filesep, 'processed', filesep, ['Sub_' sub_id{sub_counter} '_SocPriors.mat']), 'explicit');
                explicit            = sum([explitive(1).explicit, explitive(2).explicit]);
                
            end
            
        end
        
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
        
        subdata                  =  [num(subidx-1, [offer_col, acc_col, pay_col]), soc_array, opp_array, repmat(explicit, length(soc_array), 1)]; % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        
        all_sub_data(sub_counter).subdata  = subdata;
        all_sub_data(sub_counter).sub_id   = sub_id{sub_counter};

        
    end
    
end


% Go through all the data to make sure there's no duplicate subjects in the
% the multiple files that we include, and if there are take the most
% complete version for that subject 
clear all_subdata_unique
uniquesubcount = 0;
uniqueidx      = unique(sub_id);
for ii = 1:length(uniqueidx)
    repsub         = find(strcmp(sub_id, uniqueidx{ii}));
    
    for kk = 1:length(repsub)
        tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
        tmp(kk).sub_id  = uniqueidx{ii};
        all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
        all_subdata_unique(uniquesubcount).sub_id = tmp(bigidx).sub_id;
    end
    
end


for ii = 1:length(all_subdata_unique)
    sub_id          = all_subdata_unique(ii).sub_id;
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    colnames        = {'offer', 'accepted', 'payment', 'social', 'opponent', 'explicit'};
    
    flnm            = strcat('Sub_', sub_id, '_UG');
    flnm_to_save = fullfile(['..', filesep, 'data', filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name', 'colnames')
end


%% Dictator
clear all
raw_dir     = ['..', filesep, 'data' filesep,'raw'];

% social priors
pnonsoc_files    = dir([raw_dir,filesep, '*my_dictator*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(pnonsoc_files)
    txtnm = fullfile([raw_dir,filesep,pnonsoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    op1_col       = find(strcmp(txt(1, :), 'group.option1'));
    op2_col       = find(strcmp(txt(1, :), 'group.option2'));
    offer_col     = find(strcmp(txt(1, :), 'group.amount_offered'));
    pay_col       = find(strcmp(txt(1, :), 'group.payout'));
    
    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;
        %         sub_id{sub_counter}      = all_sub_name{subsub};
        tmp_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), tmp_id{sub_counter}));
        
        subdata                  = num(subidx-1, [op1_col, op2_col, offer_col, pay_col]); % subidx-1 because num doesn't read in the column headers and it therefore has 1 less row than txt
        
        subdata(any(isnan(subdata), 2), :) = []; % remove NaN rows
        sub_id{sub_counter}      = all_sub_name{subsub};
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
    
%     for kk = 1:length(repsub)
%         tmp(kk).data     = all_sub_data(repsub(kk)).subdata;
%         all_lengths(kk)  = length(all_sub_data(repsub(kk)).subdata);
%     end
%     
%     % then take the biggest one
%     if ~isempty(all_sub_data(repsub(kk)).subdata)
%         uniquesubcount = uniquesubcount +1;
%         bigidx = find(max(all_lengths));
%         all_subdata_unique(uniquesubcount).subdata  = tmp(bigidx).data;
%     else
%         fprintf('%s is an empty matrix, row %d of all_sub_data.\n', all_sub_data(repsub(kk)).sub_id, repsub(kk));
%     end
%     
%     
    for kk = 1:length(repsub)
        tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
        tmp(kk).sub_id  = uniqueidx{ii};
        all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
        all_subdata_unique(uniquesubcount).sub_id = tmp(bigidx).sub_id;
    end
    
end


for ii = 1:length(all_subdata_unique)
%     sub_id          = uniqueidx{ii};
    sub_id          = all_subdata_unique(ii).sub_id;
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('Sub_', sub_id, '_Dictator');
    flnm_to_save = fullfile(['..', filesep, 'data', filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end

%% Risk

clear all
raw_dir     = ['..', filesep, 'data' filesep,'raw'];

% social priors
pnonsoc_files    = dir([raw_dir,filesep, '*my_risk*xlsx*']);
sub_name_col  = 2;
sub_counter   = 0;

for f = 1:length(pnonsoc_files)
    txtnm = fullfile([raw_dir,filesep,pnonsoc_files(f).name]);
    
    [num,txt,raw] = xlsread(txtnm); % get the raw data file
    
    all_sub_name   = unique(txt(2:end, sub_name_col)); 
    
    risk_col       = find(strcmp(txt(1, :), 'player.risk_choice'));
    

    for subsub = 1:length(all_sub_name)
        sub_counter              = sub_counter + 1;

        tmp_id{sub_counter}      = all_sub_name{subsub};
        sub_name                 = ['00', num2str(subsub)];
        
        subidx                   = find(strcmp(txt(:, sub_name_col), tmp_id{sub_counter}));
        
        subdata                  = txt(subidx, risk_col);
        
        if ~isempty(subdata{:})
            subdata              = str2double(subdata{:}(end));
        end
        
        sub_id{sub_counter}      = all_sub_name{subsub};
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
    
%     for kk = 1:length(repsub)
%         tmp(kk).data     = all_sub_data(repsub(kk)).subdata;
%         all_lengths(kk)  = length(all_sub_data(repsub(kk)).subdata);
%     end
%     
%     % then take the biggest one
%     if ~isempty(all_sub_data(repsub(kk)).subdata)
%         uniquesubcount = uniquesubcount +1;
%         bigidx = find(max(all_lengths));
%         all_subdata_unique(uniquesubcount).subdata  = tmp(bigidx).data;
%     else
%         fprintf('%s is an empty matrix, row %d of all_sub_data.\n', all_sub_data(repsub(kk)).sub_id, repsub(kk));
%     end
    
    for kk = 1:length(repsub)
        tmp(kk).data    = all_sub_data(repsub(kk)).subdata;
        tmp(kk).sub_id  = uniqueidx{ii};
        all_lengths(kk) = length(all_sub_data(repsub(kk)).subdata);
    end
    
    % then take the biggest one
    if ~isempty(all_sub_data(repsub(kk)).subdata)
        uniquesubcount = uniquesubcount +1;
        bigidx = find(max(all_lengths));
        all_subdata_unique(uniquesubcount).subdata = tmp(bigidx).data;
        all_subdata_unique(uniquesubcount).sub_id = tmp(bigidx).sub_id;
    end
    
end


for ii = 1:length(all_subdata_unique)
%     sub_id          = uniqueidx{ii};
    sub_id          = all_subdata_unique(ii).sub_id;
    sub_name        = ['00', num2str(ii)];
    subdata         = all_subdata_unique(ii).subdata;
    
    flnm            = strcat('Sub_', sub_id, '_Risk');
    flnm_to_save = fullfile(['..', filesep, 'data', filesep, 'processed', filesep, flnm]);
    
    save(flnm_to_save,'subdata','sub_id','sub_name')
end