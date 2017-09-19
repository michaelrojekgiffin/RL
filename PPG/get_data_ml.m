clear all
close all
clc



cur_dir = pwd;
txtnm = fullfile(cur_dir,'hor_ml.xlsx');

data_dir = fullfile(cur_dir,'data_matlab');


[num,txt,raw] = xlsread(txtnm);

subjects = unique(num(:,4));


for k_sub = 1:length(subjects)
    
    i_sub = subjects(k_sub);
    
    sub_lines        = find(num(:,4)==i_sub);
    data             = num(sub_lines,:);
    role             = strcat(txt{sub_lines(1)+1,2});
    gender           = strcat(txt{sub_lines(1)+1,12});
    sex_night_before = strcat(txt{sub_lines(1)+1,26});
    birth_control    = strcat(txt{sub_lines(1)+1,27});
    food             = strcat(txt{sub_lines(1)+1,28});
    
    for k = 1:size(txt,2)
        column_ID{k}        = txt{1,k};
    end
    
    flnm = strcat('DATA_sub',num2str(i_sub));
    flnm_to_save = fullfile(data_dir,flnm);
    
    save(flnm_to_save,'data','role','gender','sex_night_before','birth_control','food','column_ID')
end



% [sub_inv, sub_rew, predprey] = sub_data_load_ml(sub_name,txtnm)