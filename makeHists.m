%
%   makeHist.m
%
% Make histograms of all features.
%
clear; clc;
close all;
%% Load Features
F = readtable("bicepstrum_f.csv");

%% Number of features and stages
stages = [
    "Sleep stage N1",
    "Sleep stage N2",
    "Sleep stage N3",
    "Sleep stage R",
    "Sleep stage W",
    ];

features = 2:11;

%% Draw the histograms

for i = 1:numel(features)
    fi = features(i);

    figure; hold on;
    for j = 1 : numel(stages)
        stage = stages(j);
        x = F{F{:,end}==stage, fi};
        histogram(x, 'Normalization','probability');
    end
    legend(stages);
    title(sprintf("%s", F.Properties.VariableNames{fi}));
end