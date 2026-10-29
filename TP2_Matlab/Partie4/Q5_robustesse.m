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
E = L+1 : N-L;
gain = @(dh) 10*log10(mean(d1(E).^2) / mean((d1(E) - dh(E)).^2));
fprintf('Apprentissage : %.1f s (homme) | %.1f s (femme)\n', length(a1)/Fs, length(a2)/Fs);

%% (a) Reference : r_d1 connu sur chaque trame
oracle = @(rx, idx) toeplitz(rx) \ autocorr_biased(d1(idx), P);
g_oracle = gain(wiener_ola(x, P, L, oracle));
fprintf('r_d1 connu (reference) : %.2f dB\n', g_oracle);

%% (b) Sensibilite : r_d1 juste en forme, faux en puissance (r_d1 -> c * r_d1)
cl = [0.25 0.5 0.8 1 1.25 2 4];
g_c = zeros(size(cl));
for m = 1:length(cl)
    g_c(m) = gain(wiener_ola(x, P, L, @(rx, idx) cl(m) * oracle(rx, idx)));
    fprintf('c = %.2f : %.2f dB\n', cl(m), g_c(m));
end

%% (c) Apprentissage : vecteurs r_d sur les trames actives du reste des enregistrements
R1 = vecteurs_rd(a1, P, L);          % homme : une colonne par trame (mat_rd1)
R2 = vecteurs_rd(a2, P, L);          % femme
fprintf('Trames d''apprentissage : %d (homme) | %d (femme)\n', size(R1,2), size(R2,2));
R1n = R1 ./ R1(1,:);                 % formes normalisees : r(k)/r(0)
R2n = R2 ./ R2(1,:);

figure;
subplot(2,1,1);
imagesc(R1n); colorbar; axis xy;
xlabel('trame'); ylabel('retard k');
title('mat\_rd1 : r_{d_1}(k)/r_{d_1}(0) pour chaque trame d''apprentissage');
subplot(2,1,2);
histogram(R1n(2,:), 40); grid on;
xlabel('r_{d_1}(1)/r_{d_1}(0)'); ylabel('nombre de trames');
title('Distribution de la correlation au retard 1');
saveas(gcf, [rep 'p4_q5_matrice.png']);

%% (d) Dictionnaires et gain en fonction du nombre de representants Nc
Ncl = [1 2 4 8 16];
g_enonce = zeros(size(Ncl));  g_mes = g_enonce;  g_mod = g_enonce;
for m = 1:length(Ncl)
    Nc = Ncl(m);
    D1 = kmeans_simple(R1,  Nc);     % methode de l'enonce : r_d1 bruts
    V  = kmeans_simple(R1n, Nc);     % formes normalisees (homme)
    U  = kmeans_simple(R2n, Nc);     % formes normalisees (femme)
    g_enonce(m) = gain(wiener_ola(x, P, L, @(rx, idx) toeplitz(rx) \ choix_enonce(rx, D1)));
    g_mes(m)    = gain(wiener_ola(x, P, L, @(rx, idx) filtre_decomp(rx, V, U, false)));
    g_mod(m)    = gain(wiener_ola(x, P, L, @(rx, idx) filtre_decomp(rx, V, U, true)));
    fprintf('Nc = %2d : enonce = %.2f dB | decomp. (R_x mesure) = %.2f dB | decomp. (R_x modele) = %.2f dB\n', ...
        Nc, g_enonce(m), g_mes(m), g_mod(m));
    if Nc == 4, D4 = D1; V4 = V; end
    if Nc == 8, V8 = V; U8 = U; end
end

figure;
semilogx(Ncl, g_enonce, 'bo-', 'LineWidth', 1.3); hold on;
semilogx(Ncl, g_mes, 'ms-', 'LineWidth', 1.3);
semilogx(Ncl, g_mod, 'r^-', 'LineWidth', 1.3);
yline(g_oracle, 'k--', 'r_{d_1} connu');
yline(0, 'k:');
grid on; xlabel('nombre de representants N_c'); ylabel('gain homme (dB)');
title(sprintf('Filtrage avec dictionnaire (P = %d, L = 40 ms, OLA)', P));
legend('argmax de la correlation (enonce)', 'decomposition, R_x mesure', ...
    'decomposition, R_x modele', 'Location', 'best');
saveas(gcf, [rep 'p4_q5_gain_Nc.png']);

figure;
semilogx(cl, g_c, 'ko-', 'LineWidth', 1.3); hold on;
semilogx(cl, -20*log10(abs(1 - cl)), 'r--');
grid on; ylim([-15 15]);
xlabel('facteur c (r_{d_1} utilise = c \times r_{d_1} vrai)'); ylabel('gain homme (dB)');
title('Sensibilite a une erreur sur la puissance de d_1');
legend('gain mesure', '-20 log_{10}|1-c|', 'Location', 'best');
saveas(gcf, [rep 'p4_q5_echelle.png']);

% Spectres des 4 representants de l'homme
figure;
nf = 1024;  f = (0:nf/2-1) * Fs / nf;
for i = 1:4
    rs = [flipud(V4(2:end,i)); V4(:,i)];        % r(-P+1) ... r(P-1)
    S  = abs(fft(rs .* bartlett(length(rs)), nf));
    plot(f, 10*log10(S(1:nf/2)), 'LineWidth', 1.2); hold on;
end
grid on; xlabel('f (Hz)'); ylabel('DSP (dB)');
title('Spectres des N_c = 4 representants de la voix d''homme');
legend('v_1', 'v_2', 'v_3', 'v_4');
saveas(gcf, [rep 'p4_q5_atomes.png']);

%% Ecoute : enonce (Nc = 4), puis decomposition avec R_x modele (Nc = 8)
sound(wiener_ola(x, P, L, @(rx, idx) toeplitz(rx) \ choix_enonce(rx, D4)), Fs);   pause(N/Fs + 0.5);
sound(wiener_ola(x, P, L, @(rx, idx) filtre_decomp(rx, V8, U8, true)), Fs);


%% ===== Fonctions locales =====
function R = vecteurs_rd(a, P, L)
% Vecteurs r_d (retards 0..P-1) des trames actives de a (pas L/2)
seuil = 0.01 * mean(a.^2);
R = [];
for s = 1:L/2:length(a)-L+1
    t = a(s:s+L-1);
    if mean(t.^2) > seuil            % trames de silence ignorees
        R(:, end+1) = autocorr_biased(t, P); %#ok<AGROW>
    end
end
end

function C = kmeans_simple(X, K)
% k-means sur les colonnes de X : retourne les K centres (colonnes de C)
rng(0);
C = X(:, randperm(size(X,2), K));
for it = 1:100
    D = sum(X.^2, 1)' + sum(C.^2, 1) - 2*(X'*C);   % distances au carre
    [~, lab] = min(D, [], 2);
    Cn = C;
    for k = 1:K
        if any(lab == k), Cn(:,k) = mean(X(:, lab == k), 2); end
    end
    if isequal(Cn, C), break; end
    C = Cn;
end
end

function rd = choix_enonce(rx, D)
% Representant le plus correle a r_x (correlation normalisee)
c = (rx' * D) ./ (norm(rx) * sqrt(sum(D.^2, 1)));
[~, i0] = max(c);
rd = D(:, i0);
end

function [rd1, rxm] = decomp(rx, V, U)
% r_x ~ a v_i + b u_j (a, b >= 0) : on garde le couple qui explique le mieux r_x
% rd1 = a v_i (part de l'homme), rxm = a v_i + b u_j (modele de r_x)
best = inf;  rd1 = zeros(size(rx));  rxm = rd1;
for i = 1:size(V,2)
    for j = 1:size(U,2)
        A = [V(:,i) U(:,j)];
        c = A \ rx;
        if any(c < 0)                % contrainte a, b >= 0 : un seul representant
            a = max(0, V(:,i)'*rx / (V(:,i)'*V(:,i)));
            b = max(0, U(:,j)'*rx / (U(:,j)'*U(:,j)));
            if norm(rx - a*V(:,i)) < norm(rx - b*U(:,j)), c = [a; 0]; else, c = [0; b]; end
        end
        res = norm(rx - A*c);
        if res < best
            best = res;
            rd1 = c(1) * V(:,i);
            rxm = A * c;
        end
    end
end
end

function w = filtre_decomp(rx, V, U, modele)
% Filtre de Wiener avec r_d1 = a v_i ; R_x mesure (modele = false) ou
% R_x = toeplitz(a v_i + b u_j) coherent avec r_d1 (modele = true)
[rd1, rxm] = decomp(rx, V, U);
if modele, w = toeplitz(rxm) \ rd1; else, w = toeplitz(rx) \ rd1; end
if any(~isfinite(w)), w = zeros(size(rx)); end   % trame vide
end
