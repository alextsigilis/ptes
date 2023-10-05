%% Examine Zero-Crossing-Rate
fprintf("Computing ZCR...");

zcr1 = [];
zcr2 = [];

for i = 1 : size(dataREM,1)
    x = cell2mat(dataREM{i,1:2});
    x = (x-mean(x))./std(x);
    y = lowpass(x,12,fs);
    zcr1 = [zcr1; zerocrossrate(y)];
end


for i = 1 : size(dataWAKE,1)
    x = cell2mat(dataWAKE{i,1:2});
    x = (x-mean(x))./std(x);
    y = lowpass(x,12,fs);
    zcr2 = [zcr2; zerocrossrate(y)];
end

figure(1); hold on;
scatter(zcr1(:,1), zcr1(:,2));
scatter(zcr2(:,1), zcr2(:,2));
title("Zero Crossing Rates of EOG channels");
xlabel("zcr of EOG-E1-M2");
ylabel("zcr of EOG-E2-M2");
legend(["REM", "WAKE"]);

fprintf("Done.\n");


%% Filter signals REM
x = cell2mat(dataREM{epoch,1:2});
x = (x-mean(x))./std(x);
y = lowpass(x, 12, fs);


% Plot result
figure(2);
sgtitle("REM Signals")
subplot(2,1,1);
plot(x);
title("Pre-filter signals");
legend(["EOGE1_M2", "EOGE2_M2"]);
subplot(2,1,2);
plot(y);
title("Filtered signals");
legend(["EOGE1_M2", "EOGE2_M2"]);

%% FFT of REM signal
X = fftshift(fft(x)); X = X(idx, :);
Y = fftshift(fft(y)); Y = Y(idx, :);

% Plot result
figure(3);
sgtitle("Powerspectrum of REM signals");
subplot(2,1,1); plot(fspan, abs(X)); legend(["EOGE1_M2", "EOGE2_M2"]); title("Pre-filter PS");
subplot(2,1,2); plot(fspan, abs(Y)); legend(["EOGE1_M2", "EOGE2_M2"]); title("Filtered PS");

% Scatter plot of Angles
figure(4);
phase = unwrap(angle(Y));
scatter(phase(:,1), phase(:,2));
title("Unwraped phases of WAKE EOG Signals");
xlabel("phase of EOG-E1-M2");
ylabel("phase of EOG-E2-M2");

%% Peaks of REM Powespectrum
P = abs(Y);
P = (P-min(P)) ./ (max(P)-min(P));

w = gausswin(64,0.5); w = w / sum(w);
P(:,1) = conv(P(:,1), w, 'same');
P(:,2) = conv(P(:,2), w, 'same');

[pks1, locs1] = findpeaks(P(:,1));
[pks2, locs2] = findpeaks(P(:,2));

% Plot result
figure(5); 
sgtitle("Peaks and Phases of REM");
subplot(2,2,1); plot(fspan, P(:,1)); hold on; plot(fspan(locs1), pks1, 'o'); title("Peaks of REM Powespectrum of EOG E1 M2"); xlabel("f (Hz)")
subplot(2,2,2);plot(fspan(locs1), unwrap(angle(Y(locs1,1))));
subplot(2,2,3); plot(fspan, P(:,2)); hold on; plot(fspan(locs2), pks2, 'o'); title("Peaks of REM Powespectrum of EOG E2 M2"); xlabel("f (Hz)")
subplot(2,2,4);plot(fspan(locs2), unwrap(angle(Y(locs2,2))));

%% ============================================================================
%  ============================================================================
%  ============================================================================
clc;

%% Filter signals WAKE
x = cell2mat(dataWAKE{epoch,1:2});
x = (x-mean(x))./std(x);
y = lowpass(x, 12, fs);

% Plot result
figure(6);
sgtitle("WAKE Signals")
subplot(2,1,1);
plot(x);
title("Pre-filter signals");
legend(["EOGE1_M2", "EOGE2_M2"]);
subplot(2,1,2);
plot(y);
title("Filtered signals");
legend(["EOGE1_M2", "EOGE2_M2"]);

%% FFT of WAKE signal
X = fftshift(fft(x)); X = X(idx, :);
Y = fftshift(fft(y)); Y = Y(idx, :);

% Plot result
figure(7);
sgtitle("Powerspectrum of WAKE signals");
subplot(2,1,1); plot(fspan, abs(X)); legend(["EOGE1_M2", "EOGE2_M2"]); title("Pre-filter PS");
subplot(2,1,2); plot(fspan, abs(Y)); legend(["EOGE1_M2", "EOGE2_M2"]); title("Filtered PS");

% Scatter plot of Angles
figure(4); hold on;
phase = unwrap(angle(Y));
scatter(phase(:,1), phase(:,2));
title("Unwraped phases of EOG Signals");
xlabel("phase of EOG-E1-M2");
ylabel("phase of EOG-E2-M2");
legend(["REM", "WAKE"])

%% Peaks of WAKE Powespectrum
P = abs(Y);
P = (P-min(P)) ./ (max(P)-min(P));

w = gausswin(64,0.5); w = w / sum(w);
P(:,1) = conv(P(:,1), w, 'same');
P(:,2) = conv(P(:,2), w, 'same');

[pks1, locs1] = findpeaks(P(:,1));

[pks2, locs2] = findpeaks(P(:,2));

% Plot result
figure(8);
sgtitle("Peaks and Phases of WAKE");
subplot(2,2,1); plot(fspan, P(:,1)); hold on; plot(fspan(locs1), pks1, 'o'); title("Peaks of REM Powespectrum of EOG E1 M2"); xlabel("f (Hz)")
subplot(2,2,2);plot(fspan(locs1), unwrap(angle(Y(locs1,1))));
subplot(2,2,3); plot(fspan, P(:,2)); hold on; plot(fspan(locs2), pks2, 'o'); title("Peaks of REM Powespectrum of EOG E2 M2"); xlabel("f (Hz)")
subplot(2,2,4);plot(fspan(locs2), unwrap(angle(Y(locs2,2))));


%% Plot scatter of phases in a signle epoch
e = 3;

x = cell2mat(dataWAKE{e, 1:2}); y = lowpass(x, 12, fs);
Y = fft(Y);
p1 = unwrap(angle(Y));

x = cell2mat(dataREM{e, 1:2}); y = lowpass(x, 12, fs);
Y = fft(Y);
p2 = unwrap(angle(Y));

figure(9);
scatter(p1(:,1), p1(:,2));
hold on;
scatter(p2(:,1), p2(:,2));
legend(["WAKE", "REM"]);


%% Plot histograms of phases on the whole dataset
fprintf("Calculating phases of sample...")

phaseWAKE = [];
for i = 1 : size(dataWAKE,1)
    x = cell2mat(dataWAKE{i,1:2});
    y = lowpass(x,12,fs);
    Y = fft(y);
    phaseWAKE = [phaseWAKE; unwrap(angle(Y))];
end

phaseREM = [];
for i = 1 : size(dataREM,1)
    x = cell2mat(dataREM{i,1:2});
    y = lowpass(x,12,fs);
    Y = fft(y);
    phaseREM = [phaseREM; unwrap(angle(Y))];
end

figure(10);
subplot(2,1,1); hold on;
histogram(phaseWAKE(:,1)); histogram(phaseREM(:,1)); legend(["WAKE", "REM"]); title("EOG-E1-M2")
subplot(2,1,2); hold on;
histogram(phaseWAKE(:,2)); histogram(phaseREM(:,2)); legend(["WAKE", "REM"]); title("EOG-E2-M2")

fprintf("Done.\n");