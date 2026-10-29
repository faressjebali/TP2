function r = autocorr_ar(s, p, P)
% Autocorrelation (retards 0 a P-1) du modele AR d'ordre p ajuste sur s
% Le modele AR : s(n) = -a(2)s(n-1) - ... - a(p+1)s(n-p) + e(n), var(e) = g
s = s(:);
if mean(s.^2) < 1e-12                % trame vide : pas de modele
    r = zeros(P,1);
    return
end
win = hamming(length(s));
[a, g] = lpc(s .* win, p);           % prediction lineaire (Yule-Walker)
g = g / mean(win.^2);                % correction de la puissance de la fenetre

nf = 1024;
S = g ./ abs(fft(a(:), nf)).^2;      % DSP du modele : g / |A(f)|^2
r = real(ifft(S));                   % autocorrelation = TF inverse de la DSP
r = r(1:P);
if any(~isfinite(r)), r = zeros(P,1); end
end
