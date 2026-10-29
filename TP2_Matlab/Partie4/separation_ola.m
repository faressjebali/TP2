function dhat = separation_ola(x, d, P, L)
% Filtre de Wiener par trames de L echantillons avec recouvrement de 50 %
% et reconstruction par overlap-add (fenetre de Hann, somme des fenetres = 1)
% Entrees : x (melange), d (source cible), P (ordre), L (longueur de trame, paire)
% Sortie  : dhat (source estimee, meme longueur que x)

x = x(:); d = d(:);
N   = length(x);
H   = L/2;                           % pas d'avancement : 50 % de recouvrement
win = hann(L, 'periodic');           % fenetres decalees de L/2 : somme = 1
dhat = zeros(N,1);

for s = 1:H:N-L+1
    idx = s : s+L-1;                 % indices de la trame
    rx = xcorr(x(idx), P-1, 'biased'); rx = rx(P:end);
    rd = xcorr(d(idx), P-1, 'biased'); rd = rd(P:end);
    w  = toeplitz(rx) \ rd;          % filtre de la trame

    deb = max(1, s-P+1);             % P-1 echantillons passes de x
    y = filter(w, 1, x(deb:idx(end)));
    y = y(end-L+1:end);
    dhat(idx) = dhat(idx) + win .* y;   % overlap-add
end
end
