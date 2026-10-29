%% Partie 2 - Question (c) : convergence de eqm(k) vers EQM_min
clear; clc; close all;

% Parametres
N     = 1e4;
n     = (1:N)';
w0    = 0.05*pi;
phi   = pi/4;
sigv2 = 0.1;

% Signal propre et signal bruite
d = sin(n*w0 + phi);
v = sqrt(sigv2)*randn(N,1);
x = d + v;

Plist  = [3 5 10];
EQM    = zeros(N, length(Plist));
EQMmin = zeros(1, length(Plist));

for m = 1:length(Plist)
    Pm  = Plist(m);
    rxm = xcorr(x, Pm-1, 'biased'); rxm = rxm(Pm:end);
    rdm = xcorr(d, Pm-1, 'biased'); rdm = rdm(Pm:end);
    wm  = toeplitz(rxm) \ rdm;
    dhm = filter(wm, 1, x);

    S = 0;
    for k = 1:N
        S        = S + (d(k) - dhm(k))^2;
        EQM(k,m) = S / k;
    end
    EQMmin(m) = rdm(1) - wm' * rdm;
    fprintf('P = %2d : EQM_min = %.4f | eqm(N) = %.4f\n', Pm, EQMmin(m), EQM(N,m));
end

% Figures 23 a 25 : une figure par P
for m = 1:length(Plist)
    figure;
    plot(1:N, EQM(:,m), 'b', 'LineWidth', 1.3); hold on;
    plot([1 N], [EQMmin(m) EQMmin(m)], 'r--', 'LineWidth', 1.3);
    set(gca, 'XScale', 'log'); grid on;
    ylim([0 3*EQMmin(m)]);
    xlabel('k'); ylabel('eqm(k)');
    title(sprintf('Convergence de eqm(k) vers EQM_{min} (P = %d)', Plist(m)));
    legend('eqm(k)', 'EQM_{min}', 'Location', 'best');
    saveas(gcf, sprintf('fig%d.png', 22+m));      % fig23, fig24, fig25
end

% Figure 26 : les trois courbes ensemble
col = lines(3);
figure; hold on;
for m = 1:length(Plist)
    h(m) = plot(1:N, EQM(:,m), 'Color', col(m,:), 'LineWidth', 1.3);
    plot([1 N], [EQMmin(m) EQMmin(m)], '--', 'Color', col(m,:), 'LineWidth', 1);
end
set(gca, 'XScale', 'log'); grid on;
ylim([0 3*max(EQMmin)]);
xlabel('k'); ylabel('eqm(k)');
title('eqm(k) pour P = 3, 5 et 10 (pointillés : EQM_{min})');
legend(h, 'P = 3', 'P = 5', 'P = 10', 'Location', 'best');
saveas(gcf, 'fig26.png');