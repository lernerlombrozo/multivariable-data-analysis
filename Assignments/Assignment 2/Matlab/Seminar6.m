%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ajay Ganesh, Applied Multivariate Data Analysis CHE 494/694, Seminar 6
% PCA SVD PCR PLS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

close all; clear all; clc;

filename = 'foodData.csv'; %Selected csv file must be in same folder as current File
FoodData = csvread(filename,1,1);
X=[FoodData(:,2) FoodData(:,4) FoodData(:,6)]; 
%X=[1 2;3 4;5 6;7 8];

X_c=[]; %centered
X_n=[]; %normalized by stdev

for i=1:size(X,2)
    xcen=X(:,i)-mean(X(:,i)); 
    xnor=xcen/std(X(:,i)); 
    X_c=[X_c xcen]; 
    X_n=[X_n xnor];
end
clear i;

% X_c=[X(:,1)-mean(X(:,1)) X(:,2)-mean(X(:,2)) X(:,3)-mean(X(:,3))];
% % Zero centered 
%
% X_n=[(X(:,1)-mean(X(:,1)))/std(X(:,1)) (X(:,2)-mean(X(:,2)))/std(X(:,2)) (X(:,3)-mean(X(:,3)))/std(X(:,3))];
% % Normalized

Com=X_c'*X_c/(length(X)-1); % Co-variance matrix 

Crm=X_n'*X_n/(length(X)-1); % Correlation matrix

[C_evtr, C_eval]=eig(Com); 
[U, S, V]=svd(X_c);

[mat_coeff,mat_score,mat_latent]=pca(X_c);

svd_principal_component_directions=V; % Same as Matlab Co-efficients 
svd_scores=X_c*V; % Same as Matlab Scores 
svd_latents=(diag(S).^2)/(length(X)-1); % Same as Matlab latents

svd_loading=V*S(1:3,:)/sqrt(length(X)-1); % The term loading is used for this qunatity 
%is because of the relation which fits in. Com=(Loading)*(Loading)' 
Com_from_loadings=svd_loading*svd_loading';

fprintf('Eigenvalues of Covariance Matrix are [%4.2f %4.2f %4.2f]\n\n',diag(C_eval)); 
fprintf('Converted singular values of X are [%4.2f %4.2f %4.2f]\n\n',(diag(S).^2)/(length(X)-1)); 
fprintf('Matlab latent values are [%4.2f %4.2f %4.2f]\n\n',diag(mat_latent));

% Compare the right singular vectors V to the eigenvectors of covaraince matrix

%%%%%%%%%%%%%% Thin SVD and Data Reconstruction %%%%%%%%%%%%% 
[Uthin, Sthin, Vthin]=svd(X_c,'econ'); % Thin SVD formulation 
RedXRank1=U(:,1)*S(1,1)*V(:,1)';
fprintf('Rank of the reduced data matrix 1 is %d\n',rank(RedXRank1)); 
RedXRank2=U(:,1:2)*S(1:2,1:2)*V(:,1:2)';
fprintf('Rank of the reduced data matrix 2 is %d\n',rank(RedXRank2)); 

%%%%%%%%%%%%%%%%%%% Biplots %%%%%%%%%%%%%%%%%%%%
biplot(mat_coeff,'scores',mat_score);
%%%%%%%%%%%%%%%%% PCR %%%%%%%%%%%%%%%%%
% Use mat_score instead of X as the regressors
%%%%%%%%%%%%%%%% PLS %%%%%%%%%%%%%%%% % Use plsregress