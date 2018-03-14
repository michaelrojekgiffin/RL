clear all
close all
clc


load Learning_All_2018_02_15
A = [squeeze(LEARN_parametersLPP(:,:,2,1)),squeeze(LEARN_parametersLPP(:,:,2,2))];

load Learning_All_2018_02_28
B = squeeze(LEARN_parametersLPP(:,:,17));

figure
for k = 1:10
    subplot(2,5,k);
    plot(A(:,k),B(:,k),'o')
end