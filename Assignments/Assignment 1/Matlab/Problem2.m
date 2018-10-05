% * Needs to be manually modiffied if you add a new variable
% **  Wont display all data if you add a new variable
clear
clc
%reading file
filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);

%Taking the desired columns
X = [Y(:,2),Y(:,4),Y(:,6)] % *
[N, m]= size(X);

%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing a humanly understandable
%distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Sample mean; mu');
mu=transpose(mean(X))

disp('Sample variance; Z');
sigma=cov(X)
[m,n]= size(sigma); %we will use this for 2d and 3d plots

disp('Sample correlation: r');
r=corr(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2D scatter plots 
figure;
t=0
for i = 1:1:m % horizontal rule
    for j = 1:1:n %vertical rule
        t=t+1;
        subplot(n,m,t);
        scatter(X(:,j),X(:,i));
        xlabel(strcat('x',num2str(i))); ylabel(strcat('x',num2str(j)));
        title('scatter plot');
    end
end
clear i j t %for workspace cleaness

%3D scatter plot **
figure;
scatter3(X(:,1), X(:,2), X(:,3));
xlabel('x1'); ylabel('x2'); zlabel('x3');
title('3D-trivariate scatter plot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B:  Expression for the confidence ellipse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms X1 X2 X3; % *
Xtest= [X1; X2; X3]; % *
transpose(Xtest-mu)*inv(sigma)*(Xtest-mu)%obtaining numerical expression for the confidence ellipse
alpha=0.1;
chi2= chi2inv(1-alpha,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B:  Axis lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D,W]= eig(sigma)
(D*chi2).^(1/2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C:  Hipothesis testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.8;
chi2= chi2inv(1-alpha,m);
F = finv(1-alpha,m,N);
x=[8; 17; 31]
N*transpose(x-mu)*inv(sigma)*(x-mu)

%Bonferroni testing:
filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);
X = [Y(:,2),Y(:,4),Y(:,6)]; % *
mu=transpose(mean(X));
cov=cov(X);
[n,p] = size(X);
var=transpose(dot(cov,eye(p)));
alpha = .05;
C2=p*(n-1)/(n-p)*finv(1-alpha,p,(n-p));
a=mu+(var*C2/n).^(1/2)
b=mu-(var*C2/n).^(1/2)


