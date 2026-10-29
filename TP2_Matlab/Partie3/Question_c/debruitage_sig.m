function [x, dhat, gain_dB, EQMmin] = debruitage_sig(d, P, sigv2)
% Debruitage d'un signal d par filtre de Wiener d'ordre P
% Entrees : d (signal propre), P (ordre), sigv2 (variance du bruit)
% Sorties : x (signal bruite), dhat (signal debruite), gain en dB, EQM_min

d = d(:);
N = length(d);

v = sqrt(sigv2)*randn(N,1);          % bruit blanc gaussien
x = d + v;                           % signal bruite

rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);
wopt = toeplitz(rx) \ rd;            % filtre de Wiener
dhat = filter(wopt, 1, x);           % signal debruite

EQMmin = rd(1) - wopt'*rd;           % EQM theorique

% Gain mesure : erreur avant (bruit) / erreur apres (residu)
Pd     = rd(1);
RSB_av = 10*log10(Pd / mean(v.^2));
RSB_ap = 10*log10(Pd / mean((d - dhat).^2));
gain_dB = RSB_ap - RSB_av;
end