%% Partie 2 - Question (f) : gain en RSB en fonction de l'ordre P
clear; clc; close all;

% Parametres
N     = 1e4;
n     = (1:N)';
w0    = 0.05*pi;
phi   = pi/4;
sigv2 = 0.1;
Plist = 5:5:50;
nP    = length(Plist);

% Signal propre et signal bruite (meme realisation pour tous les P)
d = sin(n*w0 + phi);
v = sqrt(sigv2)*randn(N,1);
x = d + v;

% Puissances du signal et du bruit (boucle for)
Sd = 0; Sv = 0;
for k = 1:N
    Sd = Sd + d(k)^2;
    Sv = Sv + v(k)^2;
end
Pd = Sd / N;
Pv = Sv / N;

% RSB avant debruitage
RSBav = 10*log10(Pd / Pv);
fprintf('Pd = %.4f | Pv = %.4f | RSB_av = %.2f dB\n', Pd, Pv, RSBav);

RSBap  = zeros(1, nP);
gain   = zeros(1, nP);
gainTh = zeros(1, nP);   % gain deduit de EQM_min, pour comparaison

for m = 1:nP
    P = Plist(m);

    % Filtre de Wiener d'ordre P
    rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
    rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
    wopt = toeplitz(rx) \ rd;
    dhat = filter(wopt, 1, x);
    EQMmin = rd(1) - wopt' * rd;

    % Puissance de l'erreur residuelle (boucle for)
    Se = 0;
    for k = 1:N
        Se = Se + (d(k) - dhat(k))^2;
    end
    Pe = Se / N;

    % RSB apres debruitage et gain
    RSBap(m)  = 10*log10(Pd / Pe);
    gain(m)   = RSBap(m) - RSBav;
    gainTh(m) = 10*log10(Pv / EQMmin);

    fprintf('P = %2d : RSB_ap = %.2f dB | gain = %.2f dB | gain (EQM_min) = %.2f dB\n', ...
        P, RSBap(m), gain(m), gainTh(m));
end

% Figure 36 : gain en fonction de P
figure;
plot(Plist, gain, 'b-o', 'LineWidth', 1.3, 'MarkerFaceColor', 'b'); hold on;
plot(Plist, gainTh, 'r--s', 'LineWidth', 1.3);
grid on; xticks(Plist);
xlabel('Ordre P'); ylabel('gain (dB)');
title('Gain en RSB en fonction de P (\sigma_v^2 = 0.1)');
legend('gain mesure', 'gain deduit de EQM_{min}', 'Location', 'best');
saveas(gcf, 'fig36.png');