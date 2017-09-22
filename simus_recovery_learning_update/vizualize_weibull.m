clear all
close all


offers  = 0:.1:10;
endow   = 10*ones(1,numel(offers));
k       = (.5:.5:3);        % predator true accepance thereshold (logit intercept)
n       = length(k);

weibullp = @(b,x) (1-exp(-((x)/b(1)).^b(2)));


for p = 1:length(k);
PA(p,:)     = weibullp([2,k(p)],offers);            % compute proba of accepting the offers given current model
EV(p,:)     = (endow - offers).* PA(p,:);             % compute EV of the offers given current model

beta1 = 2.5;
pc(p,:)      = exp(beta1.*EV(p,:)) ./ sum(exp(beta1.*EV(p,:)));   % multinomial choice function
end

figure
set(gcf,'Color',[1,1,1])
subplot(3,1,1)
plot(offers,PA)

subplot(3,1,2)
plot(offers,EV)

subplot(3,1,3)
plot(offers,pc)