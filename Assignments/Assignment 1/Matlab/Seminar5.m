clear; clc;
%reading file
filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);

X = [Y(:,2),Y(:,4),Y(:,6)];

%Y=B1+B2X1+B3X2

X=[ones(size(X(:,1))) X(:,1) X(:,2) X(:,1).*X(:,2)];
alpha=0.05;

[b,bint,r,rint,stats] = regress(y,X,alpha)