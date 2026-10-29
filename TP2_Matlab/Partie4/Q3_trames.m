clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));
rep = fullfile(dossier, 'figs_partie4', filesep);
if ~isfolder(rep), mkdir(rep); end
warning('off', 'MATLAB:nearlySingularMatrix');   % trames de silence
warning('off', 'MATLAB:singularMatrix');

[d1, d2, x, Fs] = charger_voix(8);
gain = @(d, dh) 10*log10(mean(d.^2) / mean((d - dh).^2));   % RSB_av = 0 dB

%% Gain en fonction de la longueur de trame, pour plusieurs P
Lms   = [5 10 20 40 100 200];        % longueurs de trame (ms)
Plist = [10 20 50];
G1 = nan(length(Plist), length(Lms));  G2 = G1;
for p = 1:length(Plist)
    P = Plist(p);
    for l = 1:length(Lms)
        L = round(Lms(l)*1e-3*Fs);
        if P > L/2, continue; end    % trop peu d'echantillons pour estimer P retards
        G1(p,l) = gain(d1, separation_trames(x, d1, P, L));
        G2(p,l) = gain(d2, separation_trames(x, d2, P, L));
        fprintf('P = %2d | L = %3d ms : gain homme = %.2f dB | gain femme = %.2f dB\n', ...
            P, Lms(l), G1(p,l), G2(p,l));
    end
end

% References : filtre global
[~, ~, ~, gG10] = separation_wf(x, d1, d2, 10);
[~, ~, ~, gG50] = separation_wf(x, d1, d2, 50);
fprintf('Global : P = 10 -> %.2f dB | P = 50 -> %.2f dB\n', gG10, gG50);

figure;
coul = {'b', 'm', 'k'};
for p = 1:length(Plist)
    semilogx(Lms, G1(p,:), [coul{p} 'o-'], 'LineWidth', 1.3); hold on;
end
yline(gG10, 'b--');  yline(gG50, 'k--');
grid on; xlabel('longueur de trame (ms)'); ylabel('gain homme (dB)');
title('Gain du filtrage par trames en fonction de la longueur de trame');
legend('P = 10', 'P = 20', 'P = 50', 'global P = 10', 'global P = 50', 'Location', 'best');
saveas(gcf, [rep 'p4_q3_gain_L.png']);

%% Criteres par trame (L = 20 ms, P = 20)
P = 20;  L = round(20e-3*Fs);  K = floor(length(x)/L);
d1h = separation_trames(x, d1, P, L);
RSBav = nan(K,1);  g = nan(K,1);
for i = 1:K
    idx = (i-1)*L+1 : i*L;
    p1 = mean(d1(idx).^2);  p2 = mean(d2(idx).^2);
    if p1 < 0.01*mean(d1.^2), continue; end    % homme silencieux : trame ignoree
    RSBav(i) = 10*log10(p1 / p2);
    g(i)     = 10*log10(p1 / mean((d1(idx) - d1h(idx)).^2)) - RSBav(i);
end
fprintf('L = 20 ms, P = 20 : gain moyen par trame (homme actif) = %.2f dB\n', mean(g, 'omitnan'));

t = ((1:K) - 0.5) * L / Fs;          % instant central de chaque trame (s)
figure;
subplot(2,1,1);
plot(t, RSBav, 'b.-'); grid on; ylabel('RSB avant (dB)');
title('Homme face a la femme, trame par trame (L = 20 ms, P = 20)');
subplot(2,1,2);
plot(t, g, 'r.-'); grid on; xlabel('t (s)'); ylabel('gain (dB)');
saveas(gcf, [rep 'p4_q3_par_trame.png']);

figure;
plot(RSBav, g, 'k.'); grid on;
xlabel('RSB avant de la trame (dB)'); ylabel('gain de la trame (dB)');
title('Gain par trame en fonction du RSB avant (L = 20 ms, P = 20)');
saveas(gcf, [rep 'p4_q3_gain_vs_RSB.png']);

%% Ecoute : global P = 50, puis trames (L = 20 ms, P = 20)
[d1g, ~, ~, ~] = separation_wf(x, d1, d2, 50);
N = length(x);
sound(d1g, Fs);  pause(N/Fs + 0.5);
sound(d1h, Fs);

d2h = separation_trames(x, d2, 20, round(20e-3*Fs));
pause(N/Fs + 0.5);
sound(d2h, Fs);

