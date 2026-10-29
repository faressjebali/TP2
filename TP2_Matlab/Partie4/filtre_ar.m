function [w, r1] = filtre_ar(xt, w0, p, P, K)
% Filtre de Wiener iteratif a modele AR sur une trame xt
% w0 : filtre initial ; K : nombre d'iterations ; p : ordre AR ; P : ordre du filtre
% A chaque iteration : estimees -> modeles AR -> nouveau filtre de Wiener
w = w0;  r1 = [];
for k = 1:K
    d1t = filter(w, 1, xt);          % estimee courante de l'homme
    r1  = autocorr_ar(d1t, p, P);    % modele AR de l'homme
    r2  = autocorr_ar(xt - d1t, p, P);   % modele AR de la femme (d2 = x - d1)
    w   = toeplitz(r1 + r2) \ r1;    % R_x modele, coherent avec r1
    if any(~isfinite(w)), w = w0; return; end
end
end
