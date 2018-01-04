clear all
close all
clc

offers  = 0:1:10;
Ra      = -[10 5 0];        % predator true accepance thereshold (logit intercept)
Rb      = [3,3,3];          % predator true accpetance moise (logit slope)
n       = length(Ra);

% logistic choice function
%--------------------------
logitp = @(b,x) exp(b(1)+b(2).*(x))./(1+exp(b(1)+b(2).*(x)));



figure
set(gcf,'Color',[1,1,1])
hold on
for k = 1:n
    plot(offers,logitp([Ra(k),Rb(k)],offers),'-o',...
        'Color',[k/n,0,(n-k)/n],...
        'MarkerFaceColor',[k/n,0,(n-k)/n],...
        'MarkerEdgeColor',[k/n,0,(n-k)/n])
    
    ll{k} = num2str(Ra(k));
end

legend(ll);
    
