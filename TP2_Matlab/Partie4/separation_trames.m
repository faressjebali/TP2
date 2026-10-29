function dhat = separation_trames(x, d, P, L)
% Filtre de Wiener recalcule sur chaque trame de L echantillons (sans recouvrement)
% Entrees : x (melange), d (source cible), P (ordre), L (longueur de trame)
% Sortie  : dhat (source estimee, meme longueur que x)

x = x(:); d = d(:);
N = length(x);
K = floor(N / L);                    % nombre de trames
dhat = zeros(N,1);

for i = 1:K
    idx = (i-1)*L+1 : i*L;           % indices de la trame i
    rx = xcorr(x(idx), P-1, 'biased'); rx = rx(P:end);
    rd = xcorr(d(idx), P-1, 'biased'); rd = rd(P:end);
    w  = toeplitz(rx) \ rd;          % filtre de la trame i

    % On filtre avec les P-1 echantillons passes de x (pas de remise a zero)
    deb = max(1, idx(1)-P+1);
    y = filter(w, 1, x(deb:idx(end)));
    dhat(idx) = y(end-L+1:end);
end
end
