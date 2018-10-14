% * Needs to be manually modiffied if you add a new variable
% **  Wont display all data if you add a new variable
close all; clear all; clc; 
%reading file
filename = 'dataSet.csv'; %Selected csv file must be in same folder as current File
DataSet = csvread(filename,1,0);

%Taking the desired columns
X = [DataSet(:,2),DataSet(:,3),DataSet(:,4)]; % *
Y = [DataSet(:,5),DataSet(:,6)]; % *
[sampleN, samples]= size(X);


%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing an easy distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mvregress for Multiple Multivariate Regression
load('flu')

%{
figure;
scatter3(X(:,2),X(:,3),Y(:,1),'filled')
hold on
x1fit = min(X(:,2)):0.1:max(X(:,2));
x2fit = min(X(:,3)):0.1:max(X(:,3));
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT; %Regression Model
mesh(X1FIT,X2FIT,YFIT);
xlabel('x1');
ylabel('x2');
zlabel('y1');
view(50,10);
hold off

figure;
scatter3(X(:,2),X(:,3),Y(:,2),'filled')
hold on
x1fit = min(X(:,2)):0.1:max(X(:,2));
x2fit = min(X(:,3)):0.1:max(X(:,3));
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT; %Regression Model
mesh(X1FIT,X2FIT,YFIT);
xlabel('x1');
ylabel('x2');
zlabel('y2');
view(50,10);
hold off
%}

%Confidence intervals of the estimate, default alpha = 0.05
disp('Confidence intervals of the estimate');
disp(bint);
