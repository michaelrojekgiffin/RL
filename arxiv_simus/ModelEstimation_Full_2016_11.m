% This function generates the likelihood of each model/paramters

function negLL = ModelEstimation_Full_2016_11(params,o,r,model)
% comp params
beta1 = params(1);                                                            % choice temperature
V0    = params(2);
lr1   = params(3);                                                             % supraliminal learning rate
% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ntrial = length(o);
%funct
maxp   = @(b,x) (x-median(offers)) > (b(1)-median(offers));
% pre-allocate
PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
V       = NaN(ntrial,numel(offers)); % estimated expected value of all offers
%initiate
V(1)     = V0;
lik     = 0;


negLL = 0;
for i = 1:ntrial
    
    if model==1
        PA(i,:)     = maxp(V(i,:),offers);                              % compute proba of accepting the offers given current model
        EV(i,:)     = (endow - offers).* PA(i,:);                       % compute EV of the offers given current model
                
        pc = exp(beta1.*EV(i,:)) ./ sum(exp(beta1.*EV(i,:)));          % multinomial choice function

        lik = pc(o(i)+1);
        
        V(i+1,:) = V(i,:) - lr1.*sign(r(i)-.5);                        % WARNING (if gain, youd ecrease the thereshold, if loss you increase it... hence the negative sign)
        
        negLL = negLL - log(lik);
    end
end

% lik = -lik;                                                                               % loglikelyhood vector
