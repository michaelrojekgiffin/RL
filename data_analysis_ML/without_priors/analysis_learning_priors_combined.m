clear all
close all
clc


PRIORS = load('Priors_All_2018_01_08');
LEARNING = load('Learning_All_2017_12_08.mat');
nLM = 2; % nlearning model

X = PRIORS.parametersLPP;
Y = squeeze(LEARNING.LEARN_parametersLPP([1:47 49:end],:,nLM,:));

figure
set(gcf,'Color',[1,1,1])
for kSoc = 1:2;
subplot(1,3,1)
hold on
plot(X(:,1,1),Y(:,2,1),'ob')
plot(X(:,1,2),Y(:,2,2),'or')
xlabel('prior intercept')
ylabel('learning intercept')

subplot(1,3,2)
hold on
plot(X(:,2,1),Y(:,3,1),'ob')
plot(X(:,2,2),Y(:,3,2),'or')
xlabel('prior slope')
ylabel('learning slope')


subplot(1,3,3)
hold on
plot(X(:,3,1),Y(:,1,1),'ob')
plot(X(:,3,2),Y(:,1,2),'or')
xlabel('prior temperature')
ylabel('learning temperature')
end

        
rfx = 1:100;
Soc_sub = [ones(100,1);zeros(100,1)]; % BIG vs SMALL stake
Ppr_sub = ([X(:,1,2);X(:,1,1);]);
Plearn_sub = ([Y(:,2,2);Y(:,2,1)]);
MLE_sub = repmat(rfx',2,1);
varnames = {'Soc_sub','Ppr_sub','Plearn_sub','MLE_sub'};
tbl = table(Soc_sub,Ppr_sub,Plearn_sub,MLE_sub(:),'VariableNames', varnames);
glme1 = fitglme(tbl,'Plearn_sub ~ 1 + Ppr_sub + Soc_sub + (1 |MLE_sub)','Distribution','Normal');
disp(glme1)


rfx = 1:100;
Soc_sub = [ones(100,1);zeros(100,1)]; % BIG vs SMALL stake
Ppr_sub = [X(:,2,2);X(:,2,1);];
Plearn_sub = [Y(:,3,2);Y(:,3,1)];
MLE_sub = repmat(rfx',2,1);
varnames = {'Soc_sub','Ppr_sub','Plearn_sub','MLE_sub'};
tb2 = table(Soc_sub,Ppr_sub,Plearn_sub,MLE_sub(:),'VariableNames', varnames);
glme2 = fitglme(tb2,'Plearn_sub ~ 1 + Ppr_sub + Soc_sub + (1 |MLE_sub)','Distribution','Normal');
disp(glme2)
        


rfx = 1:100;
Soc_sub = [ones(100,1);zeros(100,1)]; % BIG vs SMALL stake
Ppr_sub = [X(:,3,2);X(:,3,1);];
Plearn_sub = [Y(:,1,2);Y(:,1,1)];
MLE_sub = repmat(rfx',2,1);
varnames = {'Soc_sub','Ppr_sub','Plearn_sub','MLE_sub'};
tb3 = table(Soc_sub,Ppr_sub,Plearn_sub,MLE_sub(:),'VariableNames', varnames);
glme3 = fitglme(tb3,'Plearn_sub ~ 1 + Ppr_sub + Soc_sub + (1 |MLE_sub)','Distribution','Normal');
disp(glme3)
        
% 
% [b,stats] = robustfit(X(:,3,1),Y(:,1,1))
% 
% figure
% plot(X(:,1,2)-X(:,1,1),Y(:,2,2)-Y(:,2,1),'o')
% [Rsoc,Psoc] = corr(X(:,1,2)-X(:,1,1),Y(:,2,2)-Y(:,2,1))