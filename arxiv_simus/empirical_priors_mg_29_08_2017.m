% clear all
close all

cur_dir     = pwd;
data_dir    = fullfile(cur_dir,'data_matlab');
fl_dir      = dir(strcat(data_dir,filesep,'DATA_sub*'));
% sub_o = NaN(length(fl_dir), 1); % preallocation needs to account for it
% being just prey
k_prey = 0;
for k_sub   = 1:length(fl_dir)
    
    flnm    = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
        case 'prey'
            
            k_prey  = k_prey+1;
            sub_o(k_prey, 1)   = data(1,5);
            
    end
end

%%
sub_o_freq = zeros(1, 11);

for ii = 1:length(sub_o_freq)
    if ~isempty(sub_o(sub_o == ii-1))
        sub_o_freq(ii) = length(sub_o(sub_o == ii-1))/length(sub_o);
    else
        sub_o_freq(ii) = 0;
    end
end

