function [d1, d2, x, Fs, a1, a2] = charger_voix2(Tmax)
% Charge les deux voix (homme d1, femme d2), meme duree Tmax (s),
% centrees et de meme puissance (RSB = 0 dB), et forme le melange x.
% a1, a2 : reste des enregistrements (apres Tmax), pour l'apprentissage

dossier = fileparts(mfilename('fullpath'));
url = 'https://www.voiptroubleshooter.com/open_speech/american/';
f1  = fullfile(dossier, 'OSR_us_000_0030_8k.wav');   % voix d'homme
f2  = fullfile(dossier, 'OSR_us_000_0010_8k.wav');   % voix de femme
if ~isfile(f1), websave(f1, [url 'OSR_us_000_0030_8k.wav']); end
if ~isfile(f2), websave(f2, [url 'OSR_us_000_0010_8k.wav']); end
[y1, Fs1] = audioread(f1);
[y2, Fs2] = audioread(f2);
y1 = y1(:,1);  y2 = y2(:,1);
if Fs1 ~= Fs2
    y2 = resample(y2, Fs1, Fs2);
end
Fs = Fs1;

N  = min([length(y1), length(y2), round(Tmax*Fs)]);
d1 = y1(1:N) - mean(y1(1:N));                  % sources centrees
d2 = y2(1:N) - mean(y2(1:N));
k  = sqrt(mean(d1.^2) / mean(d2.^2));
d2 = k * d2;                                   % meme puissance -> RSB = 0 dB
x  = d1 + d2;                                  % melange

a1 = y1(N+1:end) - mean(y1(N+1:end));          % donnees d'apprentissage
a2 = k * (y2(N+1:end) - mean(y2(N+1:end)));    % meme facteur que d2
end
