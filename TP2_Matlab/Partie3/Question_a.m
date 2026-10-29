clear; clc; close all;
rng(0);
if ~isfolder('figs_partie3'), mkdir('figs_partie3'); end
%% Chargement du son
[y, Fs] = audioread('monson.wav');   % lit le son et sa frequence
y = y(:,1);                          % garde un seul canal (si stereo)

N     = 1e5;                         % nombre d'echantillons
debut = round(1*Fs);                 % saute la 1re seconde (silence)
if debut + N > length(y)             % si le son est trop court
    debut = 0;
    N = min(N, length(y));
end
d = y(debut+1 : debut+N);            % signal de reference d(n)
d = d / max(abs(d));    % amplitude ramenee entre -1 et 1
fprintf('Fs = %d Hz | N = %d | duree = %.2f s\n', Fs, N, N/Fs);

%% Question preliminaire : sigma_v^2 pour RSB = 0 dB
Pd = mean(d.^2);                     % puissance du signal
fprintf('Pd = %.4f -> sigma_v^2 pour RSB = 0 dB : %.4f\n', Pd, Pd);

%% Question (a) : EQM_min pour P = 10 et sigma_v^2 = 0.1
P     = 10;
sigv2 = 0.1;

x = d + sqrt(sigv2)*randn(N,1);      % signal bruite

% Filtre de Wiener
rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
wopt = toeplitz(rx) \ rd;            % coefficients optimaux
dhat = filter(wopt, 1, x);           % signal debruite

% EQM minimale (rd(1) correspond a r_d(0))
EQMmin = rd(1) - wopt'*rd;
fprintf('EQM_min = %.4f (bruit avant filtrage : %.2f)\n', EQMmin, sigv2);

% RSB avant / apres et gain
RSB_av = 10*log10(rd(1)/sigv2);
RSB_ap = 10*log10(rd(1)/EQMmin);
fprintf('RSB avant = %.2f dB | RSB apres = %.2f dB | gain = %.2f dB\n', ...
    RSB_av, RSB_ap, RSB_ap - RSB_av);

%% Affichage d'une fenetre
idx = 1001:1300;
figure;
plot(idx, x(idx), 'Color', [0.7 0.7 0.7]); hold on;
plot(idx, d(idx), 'b', 'LineWidth', 1.2);
plot(idx, dhat(idx), 'r', 'LineWidth', 1.2);
grid on; xlabel('n'); ylabel('Amplitude');
title(sprintf('Debruitage audio (P = %d, \\sigma_v^2 = %.1f)', P, sigv2));
legend('x(n) : signal bruite', 'd(n) : signal propre', 'dhat(n) : signal debruite');

%% Ecoute
sound(d, Fs);    pause(N/Fs + 0.5);  % son original
sound(x, Fs);    pause(N/Fs + 0.5);  % son bruite
sound(dhat, Fs);                     % son debruite

saveas(gcf, 'figs_partie3/p3_qa_signaux.png');