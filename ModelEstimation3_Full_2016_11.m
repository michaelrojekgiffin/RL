% This function generates the likelihood of each model/paramters

function negLL = ModelEstimation3_Full_2016_11(params,o,r,model, predprey)
% comp params
beta1 = params(1); % choice temperature
V0    = params(2); % initial intercept
lr1   = params(3); % supraliminal learning rate
beta2 = params(4); % slope of opponent's acceptance function?

% task param
offers  = 0:1:10;
endow   = 10*ones(1,numel(offers));
ntrial = length(o);

%funct
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));

% pre-allocate
PA      = NaN(ntrial,numel(offers)); % estimated probability of accepting all offers
EV      = NaN(ntrial,numel(offers)); % estimated expected value of all offers
% V       = NaN(ntrial); % estimated expected value of all offers
 EPc     = NaN(ntrial); % estimated expected value of all offers
 CPE     = NaN(ntrial); % estimated expected value of all offers

%initiate
V(1)     = V0;
negLL = 0;

for i = 1:ntrial
    
    if model==1
        PA(i,:)     = logitp([V(i,:),beta2],offers);            % compute proba of accepting the offers given current model
% %         EV(i,:)     = (endow - offers).* PA(i,:);               % compute EV of the offers given current model
        switch predprey
            case 'prey'
                EV(i,:)     = (endow - offers).* PA(i,:);   % compute EV of the offers given current model
            case 'predator'
                EV(i,:)     = (endow - offers) +((endow - offers).* PA(i,:));   % compute EV of the offers given current model
        end
        pc = exp(beta1.*EV(i,:)) ./ sum(exp(beta1.*EV(i,:)));   % multinomial choice function

        lik = pc(o(i)+1);

        EPc(i)      = logitp([V(i,:),beta2],o(i));
        CPE(i)      = r(i) - EPc(i);
        V(i+1,:)    = V(i,:) + lr1.*CPE((i));                   % WARNING (if gain, you decrease the thereshold, if loss you increase it... hence the negative sign)
        
        negLL = negLL - log(lik);
    end
end

% lik = -lik;                                                                               % loglikelyhood vector
