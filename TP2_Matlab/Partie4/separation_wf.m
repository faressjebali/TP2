function [dhat, EQMmin, EQMmes, gain_dB, w] = separation_wf(x, d, b, P)
% Estimation de la source d dans le melange x = d + b par filtre de Wiener
% Entrees : x (melange), d (source cible), b (autre source), P (ordre)
% Sorties : dhat (source estimee), EQM theorique et mesuree, gain en dB, w

x = x(:); d = d(:); b = b(:);

rx = xcorr(x, P-1, 'biased'); rx = rx(P:end);
rd = xcorr(d, P-1, 'biased'); rd = rd(P:end);   % r_dx = r_d (sources decorrelees)
w  = toeplitz(rx) \ rd;              % filtre de Wiener
dhat = filter(w, 1, x);              % source estimee

EQMmin = rd(1) - w'*rd;              % EQM theorique
EQMmes = mean((d - dhat).^2);        % EQM mesuree

% Gain : RSB avant (d face a b) / RSB apres (d face au residu)
RSB_av  = 10*log10(mean(d.^2) / mean(b.^2));
RSB_ap  = 10*log10(mean(d.^2) / EQMmes);
gain_dB = RSB_ap - RSB_av;
end
