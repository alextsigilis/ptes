function bic = bicoherency(x, fmax)
%BICOHERENCY Estimates the bicoherency using the Direct method.
%   bic = bicoherency(x, fmax) where `x` is a matrix with dimentions [M K]
%   and `fmax` is the maximum (normalized) frequency of `x`. Each 
%   collumn of `x`is a segment of (with M samples) and there are K segments
%   in total.
%
%   The bicoherency is defined as
%
%                       | sum X(k1) X(k2) X*(k1+k2) |
%       b(k1,k2) = -----------------------------------------
%                   [ sum |X(k1) X(k2) X(k1+k2)|^2 ]^(1/2)
%
%
[M, K] = size(x);
kmax = floor(M*fmax);
y = (x-mean(x,1)) ./ std(x,1);
Y = fft(y) .* hann(M);
Y = Y(1:kmax, :);
Y = reshape(Y, [kmax, 1, K]);
Y12 = pagemtimes(Y, pagetranspose(Y));
CY = zeros(kmax,kmax,K);
for k = 1:K
    CY(:,:,k) = hankel(Y(:,k));
end
B = mean(Y12.*CY, 3);
P = mean(abs(Y12.*CY).^2, 3).^(1/2);
bic = abs(B) ./ (P + 1e-5);
end