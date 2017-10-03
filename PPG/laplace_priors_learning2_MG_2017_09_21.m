% This function provide the computational model with paramters priors
function [post, EV, PA, V]=laplace_priors_learning2_MG_2017_09_21(params,o,r,a0,b0,nmodel, lr_upper_bound, predprey, R_o)

% log prior of parameters
switch nmodel
    case 1
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        %         lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
%         plr2    = log(betapdf(lr2,1.1,1.1));                % the parameters are distrubution with mean + variance (different shapes) (?)
        
        
    case 2
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        plr2    = log(betapdf(lr2,1.1,1.1));                % the parameters are distrubution with mean + variance (different shapes) (?)
        
        
    case 3
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        %         lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
%         plr2    = log(betapdf(lr2,1.1,1.1));                % the parameters are distrubution with mean + variance (different shapes) (?)
        
        
    case 4
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        plr2    = log(betapdf(lr2,1.1,1.1));                % the parameters are distrubution with mean + variance (different shapes) (?)
        
end

% [1 .5 .5 .5 0]
% pbeta1  = log(gampdf(beta1,1.2,5.0));
% pbeta1  = log(gampdf(beta1,1.2,8.0));
% plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
% plr2    = log(betapdf(lr2,1.1,1.1));                % the parameters are distrubution with mean + variance (different shapes) (?)

%%.5 .5 -.5
if nmodel==1
    p = [pbeta1 plr1];
elseif nmodel==2
    p = [pbeta1 plr1 plr2];
elseif nmodel==3
    p = [pbeta1 plr1];
elseif nmodel==4
    p = [pbeta1 plr1 plr2];
end

p = -sum(p);

[l, EV, PA, V] = learning_models_estim_MG_2017_09_21(params,o,r,a0,b0,nmodel, predprey, R_o);

post = p + l;


