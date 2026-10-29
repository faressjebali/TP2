clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));
rep = fullfile(dossier, 'figs_partie4', filesep);
if ~isfolder(rep), mkdir(rep); end

[d1, d2, x, Fs] = charger_voix(8);

%% Gain du filtre global en fonction de P
Plist = [2 5 10 20 30 50 75 100 150 200];
g1 = zeros(size(Plist));  g2 = zeros(size(Plist));
for m = 1:length(Plist)
    P = Plist(m);
    [~, ~, ~, g1(m)] = separation_wf(x, d1, d2, P);
    [~, ~, ~, g2(m)] = separation_wf(x, d2, d1, P);
    fprintf('P = %3d : gain homme = %.2f dB | gain femme = %.2f dB\n', P, g1(m), g2(m));
end

figure;
plot(Plist, g1, 'bo-', 'LineWidth', 1.3); hold on;
plot(Plist, g2, 'r*--', 'LineWidth', 1.1);
grid on; xlabel('P'); ylabel('gain (dB)');
title('Gain du filtre global en fonction de l''ordre P');
legend('d_1 : homme', 'd_2 : femme', 'Location', 'southeast');
saveas(gcf, [rep 'p4_q2_gain_P.png']);

%% Reponse en frequence de H1 pour plusieurs P (resolution ~ Fs/P)
Pvis = [10 50 200];
coul = {'b', 'm', 'k'};
figure;
for m = 1:length(Pvis)
    [~, ~, ~, ~, w1] = separation_wf(x, d1, d2, Pvis(m));
    [H1, fh] = freqz(w1, 1, 2048, Fs);
    plot(fh, abs(H1), coul{m}, 'LineWidth', 1.1); hold on;
end
grid on; xlabel('f (Hz)'); ylabel('|H_1(f)|');
title('Filtre extrayant l''homme pour plusieurs ordres P');
legend('P = 10', 'P = 50', 'P = 200');
saveas(gcf, [rep 'p4_q2_H1_P.png']);
