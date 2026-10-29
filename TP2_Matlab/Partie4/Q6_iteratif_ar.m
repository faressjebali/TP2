clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));
rep = fullfile(dossier, 'figs_partie4', filesep);
if ~isfolder(rep), mkdir(rep); end
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
rng(0);

[d1, d2, x, Fs, a1, a2] = charger_voix2(8);
N = length(x);
P = 50;  L = round(40e-3*Fs);        % configuration de reference
p = 12;                              % ordre AR (~ Fs/1000 + 2 a 4)
E = L+1 : N-L;
gain = @(dh) 10*log10(mean(d1(E).^2) / mean((d1(E) - dh(E)).^2));

% Dictionnaires appris sur le reste des enregistrements (section precedente)
R1 = vecteurs_rd(a1, P, L);  R1n = R1 ./ R1(1,:);
R2 = vecteurs_rd(a2, P, L);  R2n = R2 ./ R2(1,:);

%% (a) References
g_oracle = gain(wiener_ola(x, P, L, @(rx, idx) toeplitz(rx) \ autocorr_biased(d1(idx), P)));
ar_oracle = @(idx) autocorr_ar(d1(idx), p, P);
g_ar_oracle = gain(wiener_ola(x, P, L, ...
    @(rx, idx) toeplitz(ar_oracle(idx) + autocorr_ar(d2(idx), p, P)) \ ar_oracle(idx)));
fprintf('r_d1 connu : %.2f dB | modele AR de d1 et d2 connus : %.2f dB\n', g_oracle, g_ar_oracle);

%% (b) Gain en fonction du nombre d'iterations, pour deux initialisations
Kl  = 0:5;                           % K = 0 : dictionnaire seul
Ncl = [2 8];
G = zeros(length(Ncl), length(Kl));
for m = 1:length(Ncl)
    V = kmeans_simple(R1n, Ncl(m));
    U = kmeans_simple(R2n, Ncl(m));
    for k = 1:length(Kl)
        G(m,k) = gain(wiener_ola(x, P, L, ...
            @(rx, idx) filtre_ar(x(idx), filtre_decomp(rx, V, U, true), p, P, Kl(k))));
        fprintf('Nc = %d | K = %d : %.2f dB\n', Ncl(m), Kl(k), G(m,k));
    end
end

figure;
plot(Kl, G(1,:), 'bo-', 'LineWidth', 1.3); hold on;
plot(Kl, G(2,:), 'rs-', 'LineWidth', 1.3);
yline(g_oracle, 'k--', 'r_{d_1} connu');
yline(g_ar_oracle, 'm--', 'AR de d_1 et d_2 connus');
yline(4.30, 'k:', 'filtre global');
grid on; xlabel('nombre d''iterations K'); ylabel('gain homme (dB)');
title(sprintf('Wiener iteratif a modele AR (p = %d, P = %d, L = 40 ms, OLA)', p, P));
legend('initialisation N_c = 2', 'initialisation N_c = 8', 'Location', 'southeast');
saveas(gcf, [rep 'p4_q6_iterations.png']);

%% (c) Influence de l'ordre AR (Nc = 8, K = 2)
V8 = kmeans_simple(R1n, 8);  U8 = kmeans_simple(R2n, 8);
pl = [4 8 12 20 30];
g_p = zeros(size(pl));
for m = 1:length(pl)
    g_p(m) = gain(wiener_ola(x, P, L, ...
        @(rx, idx) filtre_ar(x(idx), filtre_decomp(rx, V8, U8, true), pl(m), P, 2)));
    fprintf('p = %2d (K = 2) : %.2f dB\n', pl(m), g_p(m));
end

%% (d) Exemple sur une trame ou les deux voix ont des puissances proches
H = L/2;  starts = 1:H:N-L+1;
P1 = arrayfun(@(s) mean(d1(s:s+L-1).^2), starts);
P2 = arrayfun(@(s) mean(d2(s:s+L-1).^2), starts);
ok = abs(10*log10(P1 ./ P2)) < 3;    % RSB de la trame entre -3 et +3 dB
[~, i] = max((P1 + P2) .* ok);
idx = starts(i) : starts(i)+L-1;
rx  = autocorr_biased(x(idx), P);
[w0, rd_dico] = filtre_decomp(rx, V8, U8, true);
[~, rd_ar]    = filtre_ar(x(idx), w0, p, P, 2);
fprintf('Trame exemple : t = %.2f s, RSB = %.1f dB\n', starts(i)/Fs, 10*log10(P1(i)/P2(i)));

figure;
[S, f] = spectre_r(autocorr_biased(d1(idx), P), Fs);
plot(f, S, 'k', 'LineWidth', 1.6); hold on;
plot(f, spectre_r(autocorr_biased(d2(idx), P), Fs), 'k:', 'LineWidth', 1.1);
plot(f, spectre_r(rd_dico, Fs), 'b', 'LineWidth', 1.2);
plot(f, spectre_r(rd_ar, Fs), 'r', 'LineWidth', 1.2);
grid on; xlabel('f (Hz)'); ylabel('DSP (dB)');
title(sprintf('Trame a t = %.2f s : spectre de l''homme et ses estimations', starts(i)/Fs));
legend('d_1 vrai', 'd_2 vrai', 'dictionnaire a v_i (N_c = 8)', 'AR apres K = 2 iterations', ...
    'Location', 'best');
saveas(gcf, [rep 'p4_q6_trame.png']);

%% Ecoute : dictionnaire seul (K = 0), puis K = 2 iterations (Nc = 8)
sound(wiener_ola(x, P, L, @(rx, idx) filtre_decomp(rx, V8, U8, true)), Fs);   pause(N/Fs + 0.5);
sound(wiener_ola(x, P, L, @(rx, idx) filtre_ar(x(idx), filtre_decomp(rx, V8, U8, true), p, P, 2)), Fs);


%% ===== Fonction locale =====
function [S, f] = spectre_r(r, Fs)
% DSP (dB) a partir des retards 0..P-1 d'une autocorrelation
nf = 1024;  f = (0:nf/2-1)' * Fs / nf;
rs = [flipud(r(2:end)); r(:)];       % r(-P+1) ... r(P-1)
S  = abs(fft(rs .* bartlett(length(rs)), nf));
S  = 10*log10(S(1:nf/2) + eps);
end
