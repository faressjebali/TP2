clear; clc; close all;
dossier = fileparts(mfilename('fullpath'));   % dossier du script
rep = fullfile(dossier, 'figs_partie4', filesep);
if ~isfolder(rep), mkdir(rep); end

%% Chargement des deux voix (Open Speech Repository, 8 kHz)
url = 'https://www.voiptroubleshooter.com/open_speech/american/';
f1  = fullfile(dossier, 'OSR_us_000_0030_8k.wav');   % voix d'homme
f2  = fullfile(dossier, 'OSR_us_000_0010_8k.wav');   % voix de femme
if ~isfile(f1), websave(f1, [url 'OSR_us_000_0030_8k.wav']); end  % telechargement
if ~isfile(f2), websave(f2, [url 'OSR_us_000_0010_8k.wav']); end
[d1, Fs1] = audioread(f1);
[d2, Fs2] = audioread(f2);
d1 = d1(:,1);  d2 = d2(:,1);
if Fs1 ~= Fs2                        % meme frequence d'echantillonnage
    d2 = resample(d2, Fs1, Fs2);
end
Fs = Fs1;

Tmax = 8;                            % duree gardee (s)
N  = min([length(d1), length(d2), round(Tmax*Fs)]);   % meme longueur
d1 = d1(1:N);  d2 = d2(1:N);
d1 = d1 - mean(d1);  d2 = d2 - mean(d2);       % sources centrees
d2 = d2 * sqrt(mean(d1.^2) / mean(d2.^2));     % meme puissance -> RSB = 0 dB
x  = d1 + d2;                        % melange
fprintf('Fs = %d Hz | N = %d | duree = %.2f s\n', Fs, N, N/Fs);
fprintf('Puissances : Pd1 = %.4f | Pd2 = %.4f\n', mean(d1.^2), mean(d2.^2));

%% Spectres des deux sources (base de la separation)
[S1, f] = pwelch(d1, hamming(512), 256, 1024, Fs);
S2      = pwelch(d2, hamming(512), 256, 1024, Fs);
figure;
plot(f, 10*log10(S1), 'b', 'LineWidth', 1.2); hold on;
plot(f, 10*log10(S2), 'r', 'LineWidth', 1.2);
grid on; xlabel('f (Hz)'); ylabel('DSP (dB/Hz)');
title('DSP des deux sources');
legend('d_1 : voix d''homme', 'd_2 : voix de femme');
saveas(gcf, [rep 'p4_q1_dsp.png']);

%% Filtre de Wiener global (un filtre pour tout le signal)
P = 10;
[d1h, EQM1, EQMm1, g1, w1] = separation_wf(x, d1, d2, P);
[d2h, EQM2, EQMm2, g2, w2] = separation_wf(x, d2, d1, P);
fprintf('Source 1 : EQMmin = %.4f | EQM mesuree = %.4f | gain = %.2f dB\n', EQM1, EQMm1, g1);
fprintf('Source 2 : EQMmin = %.4f | EQM mesuree = %.4f | gain = %.2f dB\n', EQM2, EQMm2, g2);

%% Reponses en frequence des deux filtres
[H1, fh] = freqz(w1, 1, 1024, Fs);
H2       = freqz(w2, 1, 1024, Fs);
figure;
plot(fh, abs(H1), 'b', 'LineWidth', 1.2); hold on;
plot(fh, abs(H2), 'r', 'LineWidth', 1.2);
plot(fh, abs(H1 + H2), 'k--');
grid on; xlabel('f (Hz)'); ylabel('|H(f)|');
title(sprintf('Reponses en frequence des filtres de Wiener (P = %d)', P));
legend('|H_1| : extrait d_1', '|H_2| : extrait d_2', '|H_1 + H_2|');
saveas(gcf, [rep 'p4_q1_freqz.png']);

%% Ecoute : melange, puis estimees
sound(x, Fs);    pause(N/Fs + 0.5);
sound(d1h, Fs);  pause(N/Fs + 0.5);
sound(d2h, Fs);
