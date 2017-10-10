clear all
close all
clc

% to be run from data_ital

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
