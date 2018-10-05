%{
sigma = [1 1 1;1 3 2;1 2 2]
chi2= chi2inv(.90,3)
[V,D,W]= eig(sigma)
(D*chi2).^(1/2)
%}

clear
clc

filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);

%Bonferroni testing:
X = [Y(:,2),Y(:,4),Y(:,6)]; % *
mu=transpose(mean(X));
cov=cov(X);
[n,p] = size(X);
var=transpose(dot(cov,eye(p)));
alpha = .05;
C2=p*(n-1)/(n-p)*finv(1-alpha,p,(n-p));
a=mu+(var*C2/n).^(1/2)