% * Needs to be manually modiffied if you add a new variable
% **  Wont display all data if you add a new variable
close all; clear; clc; 
%reading file
filename = 'dataSet.csv'; %Selected csv file must be in same folder as current File
DataSet = csvread(filename,1,0);

%Taking the desired columns
X = [DataSet(:,2),DataSet(:,3),DataSet(:,4)]; % *
Y = [DataSet(:,5),DataSet(:,6)]; % *
[n,r] = size(X);
[~,m] = size(Y);
%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing an easy distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mvregress for Multiple Multivariate Regression
X = [ones(n,1) X];
[beta,Sigma,E,CovB,logL] = mvregress(X,Y)


