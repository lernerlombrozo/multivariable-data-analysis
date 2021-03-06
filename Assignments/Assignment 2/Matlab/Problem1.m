close all; clear; clc; 
%reading file
filename = 'dataSet.csv'; %Selected csv file must be in same folder as current File
DataSet = csvread(filename,1,0);

%Taking the desired columns
x = [DataSet(:,2),DataSet(:,3),DataSet(:,4)]; % *
Y = [DataSet(:,5),DataSet(:,6)]; % *
[n,r] = size(x);
[~,m] = size(Y);
%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing an easy distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Multiple Multivariate Lineal Regression %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [ones(n,1) x];
[beta,Sigma,E,CovB,logL] = mvregress(X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B: Principal Component Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coeff,score,latent,tsquared,explained,mu] = pca(x);
xCov=cov(X);

