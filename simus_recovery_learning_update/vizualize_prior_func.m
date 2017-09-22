clear all
close all
clc

X = 0:0.1:20;
Y = gampdf(X,1.2,5.0);
Y2 = gampdf(X,1.2,8.0);
figure;
hold on
plot(X,Y)
plot(X,Y2,'r')


X = 0:0.01:3;
Y3    = (betapdf(X./3,1.1,1.1));   
figure;
plot(X,Y3)