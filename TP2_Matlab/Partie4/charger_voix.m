function [d1, d2, x, Fs] = charger_voix(Tmax)
% Charge les deux voix (homme d1, femme d2), meme duree Tmax (s),
% centrees et de meme puissance (RSB = 0 dB), et forme le melange x

dossier = fileparts(mfilename('fullpath'));
url = 'https://www.voiptroubleshooter.com/open_speech/american/';
f1  = fullfile(dossier, 'OSR_us_000_0030_8k.wav');   % voix d'homme
f2  = fullfile(dossier, 'OSR_us_000_0010_8k.wav');   % voix de femme
if ~isfile(f1), websave(f1, [url 'OSR_us_000_0030_8k.wav']); end
if ~isfile(f2), websave(f2, [url 'OSR_us_000_0010_8k.wav']); end
[d1, Fs1] = audioread(f1);
[d2, Fs2] = audioread(f2);
d1 = d1(:,1);  d2 = d2(:,1);
if Fs1 ~= Fs2
    d2 = resample(d2, Fs1, Fs2);
end
Fs = Fs1;

N  = min([length(d1), length(d2), round(Tmax*Fs)]);
d1 = d1(1:N);  d2 = d2(1:N);
d1 = d1 - mean(d1);  d2 = d2 - mean(d2);       % sources centrees
d2 = d2 * sqrt(mean(d1.^2) / mean(d2.^2));     % meme puissance -> RSB = 0 dB
x  = d1 + d2;                                  % melange
end
