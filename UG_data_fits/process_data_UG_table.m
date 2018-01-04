% to make some model free plots corroborating the statistics I've found in
% R

clear
close all force
clc

ste = @(x) (std(x))/(sqrt(length(x)));

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% find paths and import data
% this gets all subjects into one big matrix, not super elegant but closer
% to what I'm used to in R so should help speed things along
cur_dir = pwd;
%  project_name = 'RL_PreyPredator';
project_name = 'RL'; % for use in michael's dropbox
findnm = strfind(pwd,project_name);

data_dir = fullfile(cur_dir(1:findnm-1),project_name,'data','processed_R', 'UG_all.txt');
fid = fopen(data_dir);

UG_cell = textscan(fid, '%d%d%d%d%d%d%s%d%d', 'Delimiter' , '\t\n', 'Headerlines', 1);
fclose(fid);

all_sub_names = unique(UG_cell{:, 7});
sub_count = 1;
% I was thinking of adding a column to tell me opponent specific trial numbers
UG = NaN(length(UG_cell{:,1}), length(UG_cell)); 
for ii = 1:length(UG_cell{:,1})
    UG(ii, 1) = UG_cell{1,1}(ii);
    UG(ii, 2) = UG_cell{1,2}(ii);
    UG(ii, 3) = UG_cell{1,3}(ii);
    UG(ii, 4) = UG_cell{1,4}(ii);
    UG(ii, 5) = UG_cell{1,5}(ii);
    UG(ii, 6) = UG_cell{1,6}(ii);
    
    if strcmp(UG_cell{1,7}(ii), all_sub_names{sub_count})
        UG(ii, 7) = sub_count;
    else
        sub_count = sub_count + 1;
        UG(ii, 7) = sub_count;
    end
    UG(ii, 8) = UG_cell{1,8}(ii);
    UG(ii, 9) = UG_cell{1,9}(ii);
    
end

save(fullfile(cur_dir(1:findnm-1),project_name,'data','processed_R', 'UG.mat'),'UG')
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
