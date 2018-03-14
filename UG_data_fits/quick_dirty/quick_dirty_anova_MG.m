
%% FULL mixed-effect model

load('easy_sub27-Feb-2018.mat')
% (1, 1): nonsoc imp
% (1, 2): Nonsoc exp
% (2, 1): soc Imp
% (2, 2): soc exp

% for all, model parameters are in following order
% beta, a0 (intercept), b0 (slope), alpha1, alpha2

% alphadiff = [easy_sub{1,2}(:,4) - easy_sub{1,2}(:,5); easy_sub{2,2}(:,4) - easy_sub{2,2}(:,5)]; 

alphadiff(:, 1) = easy_sub{1,2}(:,4) - easy_sub{1,2}(:,5); 
alphadiff(:, 2) = easy_sub{2,2}(:,4) - easy_sub{2,2}(:,5); 


[H P CI STATS] = ttest2(alphadiff(:, 1), alphadiff(:, 2))

[H P CI STATS] = ttest2(easy_sub{2,2}(:,4), easy_sub{2,2}(:,5))
[H P CI STATS] = ttest2(easy_sub{1,2}(:,4), easy_sub{1,2}(:,5))


[H P CI STATS] = ttest2(easy_sub{2,2}(:,1), easy_sub{1,2}(:,1))

[H P CI STATS] = ttest2(easy_sub{2,2}(:,2), easy_sub{1,2}(:,2))

[H P CI STATS] = ttest2(easy_sub{2,2}(:,4), easy_sub{1,2}(:,4))

[H P CI STATS] = ttest2(easy_sub{2,2}(:,5), easy_sub{1,2}(:,5))

[H P CI STATS] = ttest2(easy_sub{2,2}(:,2), easy_sub{1,2}(:,2))

sub       = ([1:length(alphadiff)*2, 1:length(alphadiff)*2])';

alpha1     = [easy_sub{1,1}(:,4); easy_sub{1,2}(:,4);...
                easy_sub{2,1}(:,4); easy_sub{2,2}(:,4)];
            
cond    = [repmat(1, size(easy_sub{1,1}(:,4))); repmat(2, size(easy_sub{1,1}(:,4))); ...
            repmat(3, size(easy_sub{1,1}(:,4))); repmat(4, size(easy_sub{1,1}(:,4)))]; 
        
 [p_cont,table_cont,stats_cont]=anovan(alpha1,{cond,sub},'random',2,'model','full');      

%  offer = offer_mat(:);
%  sub = sub_mat(:);
%  cond = cond_mat(:);
%  soc = soc_mat(:);
%  trial = trial_mat(:);
%  vnames = {'social','cond','trial','sub'};
%  [p_cont,table_cont,stats_cont]=anovan(offer,{soc,cond,trial,sub},'random',4,'continuous',[2,3],'model','full','varnames',vnames);
% 
