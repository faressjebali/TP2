clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));
sortie  = fullfile(dossier, 'pesq_wav');
if ~isfolder(sortie), mkdir(sortie); end
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:singularMatrix');
rng(0);

[d1, d2, x, Fs, a1, a2] = charger_voix2(8);
N = length(x);
P = 50;  L = round(40e-3*Fs);  p = 12;
E = L+1 : N-L;
gain = @(d, dh) 10*log10(mean(d(E).^2) / mean((d(E) - dh(E)).^2));

R1 = vecteurs_rd(a1, P, L);  V8 = kmeans_simple(R1 ./ R1(1,:), 8);
R2 = vecteurs_rd(a2, P, L);  U8 = kmeans_simple(R2 ./ R2(1,:), 8);

%% Estimees de l'homme pour chaque methode (la femme : x - d1h)
noms = {'melange', 'global', 'trames', 'trames_ola', 'dictionnaire', 'ar_iteratif'};
D1h = cell(1, 6);
D1h{1} = x;
D1h{2} = separation_wf(x, d1, d2, P);
D1h{3} = separation_trames(x, d1, P, L);
D1h{4} = wiener_ola(x, P, L, @(rx, idx) toeplitz(rx) \ autocorr_biased(d1(idx), P));
D1h{5} = wiener_ola(x, P, L, @(rx, idx) filtre_decomp(rx, V8, U8, true));
D1h{6} = wiener_ola(x, P, L, @(rx, idx) filtre_ar(x(idx), filtre_decomp(rx, V8, U8, true), p, P, 1));

%% Ecriture des .wav (normalises pour eviter l'ecretage ; PESQ egalise les niveaux)
ecrire = @(nom, s) audiowrite(fullfile(sortie, [nom '.wav']), 0.9 * s / max(abs(s)), Fs);
ecrire('ref_homme', d1);
ecrire('ref_femme', d2);
fid = fopen(fullfile(sortie, 'gains.csv'), 'w');   % gains en RSB, lus par pesq_synthese.py
fprintf(fid, 'methode,gain_dB\n');
for m = 1:6
    d2h = x - D1h{m};
    if m == 1, d2h = x; end          % melange : aucune separation
    ecrire(['homme_' noms{m}], D1h{m});
    ecrire(['femme_' noms{m}], d2h);
    fprintf(fid, '%s,%.4f\n', noms{m}, gain(d1, D1h{m}));
    fprintf('%-13s : gain = %6.2f dB\n', noms{m}, gain(d1, D1h{m}));
end
fclose(fid);
fprintf('%d fichiers ecrits dans %s\n', numel(dir(fullfile(sortie, '*.wav'))), sortie);
fprintf('Etape suivante : python pesq_synthese.py\n');
