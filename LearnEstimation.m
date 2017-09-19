% This function generates the likelihood of each model/paramters

function negLL = LearnEstimation(params,o)
% comp params
beta0 = params(1);                                                            % choice temperature
beta1 = params(2);
beta2 = params(3);                                                            % choice temperature
beta3 = params(4);
beta4 = params(5);                                                            % choice temperature
beta5 = params(6);
betaX = params(7);

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));

%funct
logitp  = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% Compute the p(offers)
negLL = 0;
for k = 1:3
    switch k
        case 1
            PA      = logitp([beta0,beta1],offers);            % compute proba of accepting the offers given current model
        case 2
            PA      = logitp([beta2,beta3],offers);
        case 3
            PA      = logitp([beta4,beta5],offers);
    end
    
    EV      = (endow - offers).* PA;               % compute EV of the offers given current model
    pc      = exp(betaX.*EV) ./ sum(exp(betaX.*EV));
    
    %likelihood
    negLL =  negLL + sum(-log(pc(o(:,k)+1)));
end