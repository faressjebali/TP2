clear; clc; close all;
rng(0);

% Signal (comme en (a))
[y, Fs] = audioread('monson.wav');
y = y(:,1);
N = 1e5;  debut = round(1*Fs);
d = y(debut+1 : debut+N);
d = d / max(abs(d));

sigv2 = 0.1;
x = d + sqrt(sigv2)*randn(N,1);      % meme bruit pour tous les P

Plist  = 5:5:50;
EQMmin = zeros(size(Plist));
gain   = zeros(size(Plist));

% Autocorrelations calculees une seule fois jusqu'au retard max
Pmax = max(Plist);
rxF = xcorr(x, Pmax-1, 'biased'); rxF = rxF(Pmax:end);
rdF = xcorr(d, Pmax-1, 'biased'); rdF = rdF(Pmax:end);

RSB_av = 10*log10(rdF(1)/sigv2);

for m = 1:length(Plist)
    P  = Plist(m);
    rx = rxF(1:P);  rd = rdF(1:P);   % on garde les P premiers retards
    wopt = toeplitz(rx) \ rd;
    EQMmin(m) = rd(1) - wopt'*rd;
    gain(m)   = 10*log10(rd(1)/EQMmin(m)) - RSB_av;
    fprintf('P = %2d : EQM_min = %.5f | gain = %.2f dB\n', P, EQMmin(m), gain(m));
end

% (e) EQM_min en fonction de P
figure;
plot(Plist, EQMmin, 'bo-', 'LineWidth', 1.3); grid on;
xlabel('P'); ylabel('EQM_{min}');
title('EQM_{min} en fonction de P (\sigma_v^2 = 0.1)');
saveas(gcf, 'p3_qb_eqm.png');

% (f) gain en fonction de P
figure;
plot(Plist, gain, 'ro-', 'LineWidth', 1.3); grid on;
xlabel('P'); ylabel('gain (dB)');
title('Gain en RSB en fonction de P (\sigma_v^2 = 0.1)');
saveas(gcf, 'figs_partie3/p3_qb_gain.png');