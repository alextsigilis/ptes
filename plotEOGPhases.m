function plotEOGPhases(data, epoch)
x = cell2mat(data{epoch,6});
y = cell2mat(data{epoch,7});
figure;
X = fft(x);
Y = fft(y);
plot(unwrap(angle(X))); hold on; plot(unwrap(angle(Y)));
mean(abs(unwrap(angle(X))-unwrap(angle(Y))))
end