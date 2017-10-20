% This function provide the computational model with paramters priors
function [post, PE]=laplace_priors_learning2_MG(params,o,r,a0,b0,nmodel, lr_upper_bound)
switch nmodel
    case 1
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        
        %         lr2   = params(3); % supraliminal learning rate
        
    case 2
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        plr2    = log(betapdf(lr2,1.1,1.1));                                        % the parameters are distrubution with mean + variance (different shapes) (?)
        
        
    case 3
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        %         lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        
        
        
    case 4
        beta1 = params(1); % choice temperature
        lr1   = params(2); % supraliminal learning rate
        lr2   = params(3); % supraliminal learning rate
        pbeta1  = log(gampdf(beta1,1.2,8.0));
        plr1    = log(betapdf(lr1/lr_upper_bound,1.1,1.1));
        plr2    = log(betapdf(lr2,1.1,1.1));                                        % the parameters are distrubution with mean + variance (different shapes) (?)
        
end
% log prior of parameters
% beta1  = params(1);
% lr1    = params(2);                                                            % choice temperature
% lr2    = params(3);

% [1 .5 .5 .5 0]
% pbeta1  = log(gampdf(beta1,1.2,5.0));
% pbeta1  = log(gampdf(beta1,1.2,8.0));
% plr1    = log(betapdf(lr1/3,1.1,1.1));
% plr2    = log(betapdf(lr2,1.1,1.1));                                        % the parameters are distrubution with mean + variance (different shapes) (?)


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

[l PE]= learning_models_estim_MG(params,o,r,a0,b0,nmodel);

post = p + l;


