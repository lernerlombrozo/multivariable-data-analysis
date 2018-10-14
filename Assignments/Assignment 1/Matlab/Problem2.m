% * Needs to be manually modiffied if you add a new variable
% **  Wont display all data if you add a new variable
clear; clc;
%reading file
filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);

%Taking the desired columns
X = [Y(:,2),Y(:,4),Y(:,6)]; % *
[sampleN, samples]= size(X);

%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing a humanly understandable
%distribution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Sample mean; mu');
XMean=transpose(mean(X))

disp('Sample variance; Z');
XVariance=cov(X)
[m,n]= size(XVariance); %we will use this for 2d and 3d plots

disp('Sample correlation: r');
Xr=corr(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A: Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2D scatter plots 
figure;
t=0;
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
disp('numerical expression for the confidence ellipse');
sampleN*transpose(Xtest-XMean)*inv(XVariance)*(Xtest-XMean)
alpha=0.1;
disp(strcat('Chi^2(',num2str(samples),') (',num2str(alpha),')'));
chi2= chi2inv(1-alpha,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B:  Axis lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V,D,W]= eig(XVariance);
disp('eigenVectors');
W
disp('eigenValues');
D
disp('50% length across eigenVectors');
(D*chi2).^(1/2)
disp('100% length across eigenVectors');
(D*chi2).^(1/2)*2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C:  Hipothesis testing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xTest=[8; 17; 31];
alpha=0.8;
disp(strcat('F(',num2str(samples),', ',num2str(sampleN),'-',num2str(samples),')(',num2str(alpha),')'));
F = finv(1-alpha,m,sampleN)
disp('n (x-mu)^T Z^1 (x-M)');
FTest= sampleN*transpose(xTest-XMean)*inv(XVariance)*(xTest-XMean)
if FTest <= F
    disp('Reject H0 in favour of Ha');
else
    disp('Fail to reject H0 in favour of Ha');
end

%Individual testing:
var=transpose(dot(XVariance,eye(samples)));
alpha = .05;
t2=tinv((1-alpha/2),sampleN-1)
disp('Individual Minor limits');
XMean-sqrt(var/sampleN)*t2
disp('Individual Major limits');
XMean+sqrt(var/sampleN)*t2

%Bonferroni testing:
C2=samples*(sampleN-1)/(sampleN-samples)*finv(1-alpha,samples,(sampleN-samples));
disp('Bonferroni Minor limits');
XMean-(var*C2/sampleN).^(1/2)
disp('Bonferroni Major limits');
XMean+(var*C2/sampleN).^(1/2)


