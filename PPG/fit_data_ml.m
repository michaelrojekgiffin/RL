clear all
close all
clc



options     = optimset('Algorithm', 'interior-point', 'Display', 'iter-detailed', 'MaxIter', 10000);

cur_dir     = pwd;
data_dir    = fullfile(cur_dir,'data_matlab');
fl_dir      = dir(strcat(data_dir,filesep,'DATA_sub*'));


k_prey = 0;
for k_sub   = 1:length(fl_dir)
    
    flnm    = fullfile(data_dir,fl_dir(k_sub).name);
    load(flnm)
    
    switch role
        case 'predator'
        case 'prey'
            
            k_prey  = k_prey+1;
            sub_o   = data(:,5);
            sub_r   = data(:,11);
            
            n_rep           = 10;
            parameters_rep  = NaN(n_rep,3);
            ll_rep          = NaN(n_rep,1);
            
            for k_rep = 1:n_rep
                [parameters_rep(k_rep,1:3),ll_rep(k_rep,1)]=fmincon(@(x) ModelEstimation_Full_2016_11(x,sub_o,sub_r,1),[10*rand() 10*rand() rand() ],[],[],[],[],[0 0 0],[10 10 10],[],options);
            end
            
            [~,pos] = min(ll_rep);
            
            parameters(k_prey,:)    =   parameters_rep(pos(1),:);
            ll(k_prey)              =   ll_rep(pos(1),:);
            
    end
end





