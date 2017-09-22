clear all
close all
clc

logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));
logitp = @(b,x) b(3)/2+(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));
logitp = @(b,x)(1-b(3))./(1+exp(-(b(1)+b(2).*(x))));

nsamp = 10000;

offers = 0:.1:10;
endow   = 10*ones(1,numel(offers));

PA     = logitp([-10,3,.10],offers);            % compute proba of accepting the offers given current model
EV     = (endow - offers).* PA;             % compute EV of the offers given current model





beta1 = 2.5;

pc      = exp(beta1.*EV) ./ sum(exp(beta1.*EV));   % multinomial choice function


pd      = makedist('multinomial','probabilities',pc);           % estimate the pdf from the pC
kO      = random(pd,1,10000);                                           % selected offer, in the 1:numel(offer) spavce


figure
subplot(3,1,1)
plot(offers,PA)

subplot(3,1,2)
plot(offers,EV)

subplot(3,1,3)
hold on
plot(offers,pc)

kk=[];
for k = 1:numel(offers)
kk(k) = sum(kO==k)./nsamp;
end
plot(offers,kk)
