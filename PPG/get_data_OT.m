clear all
close all
clc

% to be run from data_ital

cur_dir = pwd;
txtnm = fullfile(cur_dir,'OT.xlsx');

data_dir = fullfile(cur_dir);


[num,txt,raw] = xlsread(txtnm);

subjects = unique(txt(2:end,13));


for k_sub = 1:length(subjects)
    
    i_sub = subjects(k_sub);
    
% %     sub_lines        = find(strcmp(txt(:,13), i_sub));
% %     data             = num(sub_lines,:);
% %    

    % for predator
% %     role             = strcat(txt{sub_lines(1)+1,5});
    
    role             = 'predator';
    sub_lines        = find(strcmp(txt(:,5), role) & strcmp(txt(:,13), i_sub));
    data             = num(sub_lines-1,:);
    
    flnm             = strcat('DATA_sub_',i_sub, '_', role);
    flnm_to_save     = fullfile(data_dir,flnm);
    
    save(flnm_to_save{:},'data','role')
    
    % for prey
    role             = 'prey';
    sub_lines        = find(strcmp(txt(:,5), role) & strcmp(txt(:,13), i_sub));
    data             = num(sub_lines-1,:);
    
    flnm             = strcat('DATA_sub_',i_sub, '_', role);
    flnm_to_save     = fullfile(data_dir,flnm);
    
    save(flnm_to_save{:},'data','role')
    
    
    
end
