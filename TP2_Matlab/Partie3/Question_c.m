clear; clc; close all;
rng(0);
if ~isfolder('figs_partie3'), mkdir('figs_partie3'); end

% Signal
[y, Fs] = audioread('monson.wav');
y = y(:,1);
N = 1e5;  debut = round(1*Fs);
d = y(debut+1 : debut+N);
d = d / max(abs(d));

P = 10;  sigv2 = 0.5;

% Decoupage en trames de 3 ms
L = round(3e-3 * Fs);                % 132 echantillons par trame
K = floor(N / L);                    % nombre de trames
fprintf('L = %d echantillons/trame | K = %d trames\n', L, K);

gains    = zeros(K,1);
EQMs     = zeros(K,1);
x_rec    = zeros(K*L,1);
dhat_rec = zeros(K*L,1);

% Debruitage trame par trame
for i = 1:K
    idx = (i-1)*L+1 : i*L;           % indices de la trame i
    [xi, dhi, gains(i), EQMs(i)] = debruitage_sig(d(idx), P, sigv2);
    x_rec(idx)    = xi;              % reconstruction du signal bruite
    dhat_rec(idx) = dhi;             % reconstruction du signal debruite
end

% Debruitage sur tout le signal (reference)
[x_g, dhat_g, gain_g, EQM_g] = debruitage_sig(d, P, sigv2);

fprintf('Gain global (tout le signal)  = %.2f dB\n', gain_g);
fprintf('Gain moyen par trame          = %.2f dB\n', mean(gains));

% Gain par trame
figure;
plot(1:K, gains, 'b.-'); hold on;
yline(gain_g, 'r--', 'LineWidth', 1.5);
yline(mean(gains), 'k:', 'LineWidth', 1.5);
grid on; xlabel('numero de trame i'); ylabel('gain (dB)');
title(sprintf('Gain par trame (P = %d, \\sigma_v^2 = %.1f)', P, sigv2));
legend('gain par trame', 'gain global', 'gain moyen par trame');
saveas(gcf, 'figs_partie3/p3_qc_gain_trames.png');

% Ecoute : debruitage global puis par trames
sound(dhat_g, Fs);    pause(N/Fs + 0.5);
sound(dhat_rec, Fs);