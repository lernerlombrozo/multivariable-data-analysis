%reading file
filename = 'FoodData.csv'; %FoodData.csv must be in same folder as Problem2.m
Y = csvread(filename,1,1);

%Taking the desired columns
X = [Y(:,2),Y(:,4),Y(:,6)]

%X = [ 1 2 -1; 2 4 -2; 3 6 -3] %Uncomment for testing a humanly understandable
%distribution

disp('Sample mean; mu');
mu=mean(X)

disp('Sample variance; Z');
sigma=cov(X)
[m,n]= size(sigma); %we will use this for 2d and 3d plots

disp('Sample correlation: r');
r=corr(X)

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

%3D scatter plots
%{

figure;
t=0
for i = 1:1:m % horizontal rule
    for j = 1:1:n %vertical rule
        t=t+1;
        subplot(n,m,t);
        
    end
end
clear i j t %for workspace cleaness

%}

%Hypothesis testing; obtaining a confiddence elipse
syms X1 X2 X3; %Naming variables (Not programing amount of variables right now)
Xtest= [X1; X2; X3];
transpose(Xtest-mu)*inv(sigma)*(Xtest-mu)%obtaining numerical expression for the confidence ellipse
alpha=0.1;
chi2inv(1-alpha,m);
