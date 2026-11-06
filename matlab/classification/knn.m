%
%   knn.m
%
% This script loads the dataset of features and trains a simple kNN model
% for sleep stage classification.
%
%% clear workspace
clear;
close all;
clc;

%% Load Data
fprintf("Loading data...");
file = uigetfile("*.csv");
data = readtable(file);
fprintf("Done.\n");

%% Remove missing values
mis_idx = ismissing(data);
data = data(~any(mis_idx, 2), :);

%% Perform PCA
fprintf("Doing PCA...");
X = data{:,3:10};

X = (X-mean(X)) ./ std(X);

[coeff,score,latent] = pca(X);
figure; plot(latent, 'o-');


Y = X * coeff(:,1:4);
Y = (Y-mean(Y)) ./ std(Y);

data = [array2table(Y) data(:,end)];


fprintf("Done.\n");

%% Plot class frequencies
tbl = tabulate(data.annotation);
t = cell2table(tbl,'VariableNames', ...
    {'Value','Count','Percent'});
t.Value = categorical(t.Value);
figure;
bar(t.Value,t.Count)
xlabel('Sleep stage')
ylabel('Number of occurances')

%% Split Datasets
% Select only relevant collumns
data = data(:, 3:end);
% Number of obervations
N = size(data,1);
% test size to N ratio
test_size_ratio = 0.3;
% split
test_idx = rand(N,1) <= test_size_ratio;
dsTrain = data(~test_idx, :);
dsTest = data(test_idx, :);

%% Train Classifier
fprintf("Training classifier...");

cost = 1 - eye(5);
cost(:,2) = 1.27;
cost(2,2) = 0;
cost(4,[2,3]) = 1.5;

Mdl = fitcknn(dsTrain, 'annotation', 'NumNeighbors', 16, 'Distance','euclidean', 'Cost', cost);
fprintf("Done.\n");

%% Plot results
fprintf("Creating confusion charts...");

figure;
labels = predict(Mdl, dsTrain);
confusionchart(dsTrain.annotation, labels, 'Normalization', 'row-normalized');
title("Confusion Matrix for the Training Dataset");

figure;
labels = predict(Mdl, dsTest);
confusionchart(dsTest.annotation, labels, 'Normalization', 'row-normalized');
title("Confusion Matrix for the Testing Dataset");

fprintf("Done.\n");