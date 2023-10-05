%
%
clear; clc;

%% Get the dataset dir
dir = uigetdir(".");
cd(dir);

%% Load a file and get the EOG signals
fprintf("Loading EDFs...\n");
pids = randsample(150, 10);
data = [];
for i = 1 : 10
    fprintf("Loading Patient %d...", pids(i));
    data = [data; loadEDF(pids(i))];
    fprintf("Done.\n")
end
data = data(:, ["EOGE1_M2", "EOGE2_M2", "Annotations"]);

fprintf("Loading Complete.\n");

%% Take sample dataset
rem_stage_idx = data.Annotations == "Sleep stage R";
wake_stage_idx = data.Annotations == "Sleep stage W";
nrem_stage_idx = (~rem_stage_idx) & (~wake_stage_idx);

dataREM = datasample(data(rem_stage_idx,:), 100);
dataNREM = datasample(data(nrem_stage_idx, :), 100);
dataWAKE = datasample(data(wake_stage_idx, :), 100);

%% Define constants
epoch = 4;
fs = 256;
n = 7680;
fspan = fs/n*(-n/2:n/2-1);
idx = (fspan > 2.5) & (fspan <= 49);
fspan = fspan(idx);

%% Compute Phases of wake stage
fprintf("Computing phases for WAKE stage...");
phaseW = [];
for i = 1 : size(dataWAKE,1)
    x = cell2mat(dataWAKE{i,1:2});
    x = lowpass(x, 12, fs);
    x = (x-mean(x))./std(x);
    X = fft(x);
    phaseW = [phaseW; unwrap(angle(X))];
end

fprintf("Done.\n");

%% Compute Phases of non REM stage
fprintf("Computing phases for NON-REM stage...");
phaseN = [];
for i = 1 : size(dataNREM,1)
    x = cell2mat(dataNREM{i,1:2});
    x = lowpass(x, 12, fs);
    x = (x-mean(x))./std(x);
    X = fft(x);
    phaseN = [phaseN; unwrap(angle(X))];
end
fprintf("Done.\n");

%% Compute Phases of non REM stage
fprintf("Computing phases for REM stage...");
phaseR = [];
for i = 1 : size(dataREM,1)
    x = cell2mat(dataREM{i,1:2});
    x = lowpass(x, 12, fs);
    x = (x-mean(x))./std(x);
    X = fft(x);
    phaseR = [phaseR; unwrap(angle(X))];
end
fprintf("Done.\n");

