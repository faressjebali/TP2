%% Partie 2 - Question (d) : influence de sigv2 pour P = 10
clear; clc; close all;
if ~isfolder('figs_partie2'), mkdir('figs_partie2'); end

% Parametres
N     = 1e4;
n     = (1:N)';
w0    = 0.05*pi;
phi   = pi/4;
P     = 10;
sigList = [0.01 0.05 0.1 0.5 1];
nS    = length(sigList);

% Signal propre
d = sin(n*w0 + phi);

EQM    = zeros(N, nS);
EQMmin = zeros(1, nS);

for s = 1:nS
    % Signal bruite pour cette valeur de sigv2
    x = d + sqrt(sigList(s))*randn(N,1);

    % Filtre de Wiener
    rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
    rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
    wopt = toeplitz(rx) \ rd;
    dhat = filter(wopt, 1, x);

    % eqm(k) calculee en ligne
    S = 0;
    for k = 1:N
        S        = S + (d(k) - dhat(k))^2;
        EQM(k,s) = S / k;
    end
    EQMmin(s) = rd(1) - wopt' * rd;
    fprintf('sigv2 = %.2f : EQM_min = %.4f | eqm(N) = %.4f\n', sigList(s), EQMmin(s), EQM(N,s));
end

% Figures 27 a 31 : une figure par sigv2
for s = 1:nS
    figure;
    plot(1:N, EQM(:,s), 'b', 'LineWidth', 1.3); hold on;
    plot([1 N], [EQMmin(s) EQMmin(s)], 'r--', 'LineWidth', 1.3);
    set(gca, 'XScale', 'log'); grid on;
    ylim([0 3*EQMmin(s)]);
    xlabel('k'); ylabel('eqm(k)');
    title(sprintf('eqm(k) pour P = 10 et \\sigma_v^2 = %.2f', sigList(s)));
    legend('eqm(k)', 'EQM_{min}', 'Location', 'best');
    saveas(gcf, sprintf('figs_partie2/fig%d.png', 26+s));      % fig27 a fig31
end

% Figure 32 : comparaison (axes log)
col = lines(nS);
figure; hold on;
for s = 1:nS
    h(s) = plot(1:N, EQM(:,s), 'Color', col(s,:), 'LineWidth', 1.3);
    plot([1 N], [EQMmin(s) EQMmin(s)], '--', 'Color', col(s,:), 'LineWidth', 1);
end
set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
ylim([0.5*min(EQMmin) 3*max(EQMmin)]);
xlabel('k'); ylabel('eqm(k)');
title('eqm(k) pour P = 10 et différentes valeurs de \sigma_v^2');
legend(h, arrayfun(@(v) sprintf('\\sigma_v^2 = %.2f', v), sigList, 'UniformOutput', false), 'Location', 'best');
saveas(gcf, 'figs_partie2/fig32.png');