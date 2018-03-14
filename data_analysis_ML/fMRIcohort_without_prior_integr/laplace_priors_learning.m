% This function provide the computational model with paramters priors
function [post]=laplace_priors_learning(params,o,r,nmodel)

% log prior of parameters
beta1   = params(1);    pbeta1  = log(gampdf(beta1,1.2,8.0));
a0      = params(2);    pa0     = log(normpdf(a0,-3,2));
b0      = params(3);    pb0     = log(gampdf(b0,1.2,8.0));
lr1     = params(4);    plr1    = log(betapdf(lr1,1.1,1.1));
lr2     = params(5);    plr2    = log(betapdf(lr2,1.1,1.1));

% the parameters are distrubution with mean + variance (different shapes) (?)
switch nmodel
    case 1
        p = [pbeta1 pa0 pb0 plr1 plr2];
    case 2
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S];
    case 3
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S];
    case 4
        b0S      = params(6);    pb0S     = log(gampdf(b0S,1.2,8.0));
        p = [pbeta1 pa0 pb0 plr1 plr2 pb0S];
    case 5
        lr1S     = params(6);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 plr1S];
    case 6
        lr2S     = params(6);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 plr2S];
    case 7
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S];
    case 8
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pb0S];
    case 9
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S plr1S];
    case 10
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        lr2S     = params(7);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S plr2S];
    case 11
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S pb0S];
    case 12
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S plr1S];
    case 13
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        lr2S     = params(7);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S plr2S];
    case 14
        b0S      = params(6);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pb0S plr1S];
    case 15
        b0S      = params(6);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr2S     = params(7);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pb0S plr2S];
    case 16
        lr1S     = params(6);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(7);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 plr1S plr2S];
    case 17
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(8);    pb0S     = log(gampdf(b0S,1.2,8.0));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S pb0S];
    case 18
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S plr1S];
    case 19
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S plr2S];
    case 20
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pb0S plr1S];
    case 21
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pb0S plr2S];
    case 22
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S plr1S plr2S];
    case 23
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S pb0S plr1S];
    case 24
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S pb0S plr2S];
    case 25
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S plr1S plr2S];
    case 26
        b0S      = params(6);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(7);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(8);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pb0S plr1S plr2S];
    case 27
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(8);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(9);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S pb0S plr1S];
    case 28
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(8);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr2S     = params(9);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S pb0S plr2S];
    case 29
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(9);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S plr1S plr2S];
    case 30
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(9);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pb0S plr1S plr2S];
    case 31
        a0S      = params(6);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(7);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(8);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(9);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pa0S pb0S plr1S plr2S];
    case 32
        beta1S   = params(6);    pbeta1S  = log(gampdf(beta1S,1.2,8.0));
        a0S      = params(7);    pa0S     = log(normpdf(a0S,-3,2));
        b0S      = params(8);    pb0S     = log(gampdf(b0S,1.2,8.0));
        lr1S     = params(9);    plr1S    = log(betapdf(lr1S,1.1,1.1));
        lr2S     = params(10);    plr2S    = log(betapdf(lr2S,1.1,1.1));
        p = [pbeta1 pa0 pb0 plr1 plr2 pbeta1S pa0S pb0S plr1S plr2S];
end

p = -sum(p);

l = learning_models_estim(params,o,r,nmodel);

post = p + l;


