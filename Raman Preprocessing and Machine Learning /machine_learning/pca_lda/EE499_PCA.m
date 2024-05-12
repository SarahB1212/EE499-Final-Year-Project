%% EE499 - Machine Learning PCA and LDA
% Loading data
clear all 
close all
file_path1 = 'PMMA_matfiles/PMMA_Calibrated_Dataset.mat';
file_path2 = 'PVC_matfiles/PVC_Calibrated_Dataset.mat';
file_path3 = 'PS_matfiles/PS_Calibrated_Dataset.mat';
files = {file_path1 file_path2 file_path3};
number_of_classes = 3;
classes = {'PMMA', 'PVC', 'PS'};

number_of_components = 5; % Number of principal components to keep

data = [];
for i = 1:3
    file_path = load(files{i});
    spectra = file_path.dataset;
    wn = file_path.wavenumber_axis;
    data = [data spectra];
   
end

%% Step 1 - Create Datasets
% Defining colour labels for each type of plastic
% PMMA = 199 spectra
% PVC = 199 spectra
% PS = 201 spectra
colour_labels = [ones(199,1) ;2*ones(199,1) ;3*ones(201,1)];

% Each row is a spectrum and each column represents a variable (i.e intensity)
data = data';    

save('PCA_Dataset','data', 'colour_labels');

train_percentage = 0.8;  % 80% for training, 20% for testing

% A cross-validation partition for splitting the data
seed = 100;
rng(seed);
cv = cvpartition(colour_labels, 'HoldOut', 1 - train_percentage);

% Indices for training and testing data
train_idx = cv.training;
test_idx = cv.test;

% Split the dataset and labels into training and testing sets
X_train = data(train_idx, :);
Y_train = colour_labels(train_idx);

X_test = data(test_idx, :);
Y_test = colour_labels(test_idx);


%% Step 2 - Perform PCA on Training Data
% Normalise data before performing PCA
featureMatrix = normalize(X_train); 

[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(featureMatrix);

% The scatter plot of the first 3 principal components
figure;
scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),[],Y_train,'filled');
xlabel('PC1','Fontsize',14);
ylabel('PC2','Fontsize',14);
zlabel('PC3','Fontsize',14);
%title('Principal Component Analysis of Microplastic Raman Spectra');
title('Projections onto the first 3 principal components','Fontsize',16);
colormap(jet(max(colour_labels))); % Adj st colourmap based on the maximum label value
colour_bar = colorbar('Ticks',[]); % Display colour bar
annotation('textbox', [0.88, 0.15, 0.1, 0.1], 'String', 'PMMA','EdgeColor', 'none');
annotation('textbox', [0.88, 0.45, 0.1, 0.1], 'String', 'PVC','EdgeColor', 'none');
annotation('textbox', [0.88, 0.70, 0.1, 0.1], 'String', 'PS','EdgeColor', 'none');

% The scatter plot of the first 2 principal components
figure;
scatter(SCORE(:,1),SCORE(:,2),[],Y_train,'filled');
grid on;
xlabel('PC1','Fontsize',14);
ylabel('PC2','Fontsize',14);

% The plot of the variance explained 
figure;
hold on
grid on
[charts, ax] = pareto(EXPLAINED);
charts(2).Color = [1 0 0];
charts(2).LineWidth = 1.5;
charts(2).Marker = '.';
charts(2).MarkerSize = 15;

xdata = charts(1).XEndPoints;
ydata = charts(1).YEndPoints;
for i = 1:8
    string = sprintf('%.2f %%',ydata(1,i));
    text(xdata(1,i)-0.2,ydata(1,i)+1.8,string);
end

title('Explained Variance by Different Principal Components')
xlabel('Principal Component');
ylabel('Explained Variance');
ylim([0,115])
legend('Individual Explained Variance','Cumulative Explained Variance','Location','northwest')
set(gca,'Fontsize',14);

%%
% The plot of the first 2 principal component loadings
figure;
hold on;
plot(wn,COEFF(:,1)+0.15, LineWidth=1.5);
plot(wn,COEFF(:,2),LineWidth=1.5);
plot(wn, ones(1024,1)*0.15,Color='k');
plot(wn, ones(1024,1)*0,Color='k');
set(gca,'XLim',[200 1800]);
set(gca,'YTick',[]);
xlabel('Wavenumber (cm^-^1)');
ylabel('PCA Coefficient');
set(gca,'Fontsize',14);
legend('PC1','PC2')


%% Step 3 - Perform LDA
% Using first 3 princpal components to train the LDA model
X_PCA = SCORE(:, 1:number_of_components);
LDA_Model = fitcdiscr(X_PCA, Y_train);
save('LDA_Model.mat', 'LDA_Model');

%% Step 4 - Prediction using LDA model
% Normalise test data before performing PCA
testFeatureMatrix = normalize(X_test); 
[COEFF, SCORE, LATENT,TSQUARED, EXPLAINED] = pca(testFeatureMatrix);

% Predict class labels for the test data using the first 3 princpal components
PCA_data = SCORE(:, 1:number_of_components);
predictedLabels = predict(LDA_Model, PCA_data);

% Generating confusion matrix
[CM, order] = confusionmat(Y_test,predictedLabels);

% Transposing to Predicted Classes in Rows and True Classes in Columns
CM = CM'

for i = 1:number_of_classes
    % For each class calculate the sensitivity and specificity 
    True_Postive = CM(i,i);
    False_Negative = sum(sum(CM([1:i-1, i+1:end],i)));
    True_Negative = sum(CM,'all') - sum(sum(CM(:,i))) - sum(sum(CM(i,:))) + CM(i,i);
    False_Positive = sum(sum(CM(i,:))) - CM(i,i);
    
    accuracy(i) = (True_Postive + True_Negative)/(True_Postive + True_Negative + False_Negative + False_Positive);
    sensitivity(i) = True_Postive/(True_Postive + False_Negative);
    specificity(i) = True_Negative/(True_Negative + False_Positive);
    precision(i) = True_Postive/(True_Postive + False_Positive);
    F1_score(i) = 2.*((sensitivity(i).*precision(i))/(sensitivity(i)+precision(i)));
end

fprintf('Class\tAccuracy\tSensitivity\tSpecificity\tPrecision\tF1-Score \n');
for i = 1:number_of_classes
    fprintf('%s\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n',classes{i},accuracy(i),sensitivity(i),specificity(i),precision(i),F1_score(i));
end

%% Step 5 - Prediction using LDA model with indivially calibrated spectra
clear all
close all
clc
load LDA_Model 
number_of_components = 5;
number_of_classes = 3;
classes = {'PMMA', 'PVC', 'PS'};

data = load('ML_Data/ML_Data_60_Processed.mat');
dataset = data.dataset;
dataset = dataset';
% Defining labels for each type of plastic
% PMMA = 20 spectra
% PVC = 20 spectra
% PS = 20 spectra
trueLabels = [ones(20,1) ;2*ones(20,1) ;3*ones(20,1)];

testFeatureMatrix = normalize(dataset); 

[COEFF, SCORE, LATENT,TSQUARED, EXPLAINED] = pca(testFeatureMatrix);

% Predict class labels for the test data using the first 3 princpal components
PCA_data = SCORE(:, 1:number_of_components);
predictedLabels = predict(LDA_Model, PCA_data);

% Generating confusion matrix
[CM, order] = confusionmat(trueLabels,predictedLabels);
% Transposing to Predicted Classes in Rows and True Classes in Columns
CM = CM'

for i = 1:number_of_classes
    % For each class calculate the sensitivity and specificity 
    True_Postive = CM(i,i);
    False_Negative = sum(sum(CM([1:i-1, i+1:end],i)));
    True_Negative = sum(CM,'all') - sum(sum(CM(:,i))) - sum(sum(CM(i,:))) + CM(i,i);
    False_Positive = sum(sum(CM(i,:))) - CM(i,i);
    
    accuracy(i) = (True_Postive + True_Negative)/(True_Postive + True_Negative + False_Negative + False_Positive);
    sensitivity(i) = True_Postive/(True_Postive + False_Negative);
    specificity(i) = True_Negative/(True_Negative + False_Positive);
    precision(i) = True_Postive/(True_Postive + False_Positive);
    F1_score(i) = 2.*((sensitivity(i).*precision(i))/(sensitivity(i)+precision(i)));
end

fprintf('Class\tAccuracy\tSensitivity\tSpecificity\tPrecision\tF1-Score \n');
for i = 1:3
    fprintf('%s\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n',classes{i},accuracy(i),sensitivity(i),specificity(i),precision(i),F1_score(i));
end
