%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajay Ganesh, Applied Multivariate Data Analysis CHE 494/694, Seminar 3
% Bivariate data analysis with confidence ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;
clc;

filename = 'FoodData.csv';
X = csvread(filename,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Sample mean');
X_mean=mean(X)

disp('Sample variance');
X_cov=cov(X)

disp('Sample correlation');
X_corr=corr(X)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
scatter3(X(:,2),X(:,4),X(:,6));
xlabel('x1'); ylabel('x2'); zlabel('x3');
title('3D-trivariate scatter plot');

figure;
subplot(3,3,1);
scatter(X(:,1),X(:,1));
xlabel('x1'); ylabel('x1');
subplot(3,3,2);
scatter(X(:,2),X(:,1));
xlabel('x2'); ylabel('x1');
subplot(3,3,3);
scatter(X(:,3),X(:,1));
xlabel('x3'); ylabel('x1');
subplot(3,3,4);
scatter(X(:,1),X(:,2));
xlabel('x1'); ylabel('x2');
subplot(3,3,5);
scatter(X(:,2),X(:,2));
xlabel('x2'); ylabel('x2');
subplot(3,3,6);
scatter(X(:,3),X(:,2));
xlabel('x3'); ylabel('x2');
subplot(3,3,7);
scatter(X(:,1),X(:,3));
xlabel('x1'); ylabel('x3');
subplot(3,3,8);
scatter(X(:,2),X(:,3));
xlabel('x2'); ylabel('x3');
subplot(3,3,9);
scatter(X(:,3),X(:,3));
xlabel('x3'); ylabel('x3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Bivariate Analaysis
%%Fitting a bivariate gaussian PDF to the data
dof=2; %p=2
alpha1=0.05; %95% confidence
alpha2=0.5; %50% confidence.

X_pair=[X(:,2) X(:,4)];

figure;
scatter(X_pair(:,1),X_pair(:,2))
xlabel('x1'); ylabel('x2');
title('Bivariate scatter plot');

X_pair_mean=mean(X_pair)
X_pair_cov=cov(X_pair)

k=0.5;
x1 = -k*X_pair_cov(1,1):.2:k*X_pair_cov(1,1); %values from -0.5COV(1,1) to 0.5COV(1,1) in a frequency of 0.2
x2 = -k*X_pair_cov(2,2):.2:k*X_pair_cov(2,2);

x1=x1+X_pair_mean(1);
x2=x2+X_pair_mean(2);

[X1,X2] = meshgrid(x1,x2);
F_x1x2 = mvnpdf([X1(:) X2(:)],X_pair_mean, X_pair_cov);
F_x1x2 = reshape(F_x1x2,length(x2),length(x1));

figure;
surf(x1,x2,F_x1x2);
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');
title(sprintf('Bivariate distribution of [x1 x2] with mu=[%4.2f %4.2f] and sigma=[%4.2f %4.2f;%4.2f %4.2f]',X_pair_mean(1),X_pair_mean(2),X_pair_cov(1,1),X_pair_cov(1,2),X_pair_cov(2,1),X_pair_cov(2,2)));%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% Compute the contour value from the chi2 distribution %%%%
%%% f(x1,x2)=exp(-0.5*chi2_p(alpha))/(2*pi*det(Sigma))
contour_value_0_05=exp(-0.5*chi2inv(1-alpha1,dof))/(2*pi*sqrt(det(X_pair_cov)));
contour_value_0_5=exp(-0.5*chi2inv(1-alpha2,dof))/(2*pi*sqrt(det(X_pair_cov)));

figure;
scatter(X_pair(:,1),X_pair(:,2));
xlim([x1(1) x1(end)]);
ylim([x2(1) x2(end)]);
hold on;
contour(x1,x2,F_x1x2,[contour_value_0_05 contour_value_0_05],'color','r');
contour(x1,x2,F_x1x2,[contour_value_0_5 contour_value_0_5],'color','b');
legend('Data Points','\alpha=0.05 (95% confidence)','\alpha=0.5 (50% confidence)','Location','northoutside','Orientation','vertical')
xlabel('x1');ylabel('x2');
title(sprintf('Confidence ellipse for alpha=%4.2f and alpha=%4.2f',0.05,0.5));
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%