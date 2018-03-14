clear all
close all
clc


load Learning_All_2017_12_08
A = [squeeze(LEARN_parametersLPP([1:47 49:101],:,2,1)),squeeze(LEARN_parametersLPP([1:47 49:101],:,2,2))];

load Learning_All_2018_02_27
B = squeeze(LEARN_parametersLPP(:,:,7));

figure
for k = 1:10
    subplot(2,5,k);
    plot(A(:,k),B(:,k),'o')
end