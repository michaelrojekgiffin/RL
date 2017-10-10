clear all
close all
clc



cur_dir = pwd;
txtnm = fullfile(cur_dir,'old_hor.xlsx');

data_dir = fullfile(cur_dir);


[num,txt,raw] = xlsread(txtnm);

subjects = unique(num(:,12));


for k_sub = 1:length(subjects)
    
    i_sub = subjects(k_sub);
    
    % for this dataset each subject plays as both predator and prey in
    % blocks of 60, so I'm saving two files for them, one as predator and
    % one as prey
    sub_lines        = find(num(:,12)==i_sub);
    data             = num(sub_lines(1:60),:);
    role             = strcat(txt{sub_lines(1)+1,2});
    
    for k = 1:size(txt,2)
        column_ID{k}        = txt{1,k};
    end
    
    flnm = strcat('DATA_sub',num2str(i_sub),'_',role);
    flnm_to_save = fullfile(data_dir,flnm);
    
    save(flnm_to_save,'data','role','column_ID')
    
    % and here the other role
    sub_lines        = find(num(:,12)==i_sub);
    data             = num(sub_lines(1:60),:);
    role             = strcat(txt{sub_lines(61)+1,2});
    
    for k = 1:size(txt,2)
        column_ID{k}        = txt{1,k};
    end
    
    flnm = strcat('DATA_sub',num2str(i_sub),'_',role);
    flnm_to_save = fullfile(data_dir,flnm);
    
    save(flnm_to_save,'data','role','column_ID')
end
