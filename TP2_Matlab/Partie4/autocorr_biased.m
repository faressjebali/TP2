function r = autocorr_biased(s, P)
% Autocorrelation biaisee de s pour les retards 0 a P-1 (vecteur colonne)
r = xcorr(s(:), P-1, 'biased');
r = r(P:end);
end
